// Combined VTK & Voro++ source codes
// 2-11-17 - mess created, can just fill a pre-defined shape (sphere, cylinder, etc)
// 2-12-17 - pass in STL, voro will fill inside that
// 3-8-17  - added 2nd slightly expanded shell for better point-inside detection (later removed)
// 4-9-17  - adding avg.neighbor inner.normal distance calc for cell center rejection (TODO: add tubesfromsplines, pipe diameter)
//

// flow is:
//  1. user passes in input STL, and packing-cube ("cell centers" x,y,z file)
//  2. constrcut voro wall using PointInsideSTL(x,y,z) & CutCell(x,y,z) callbacks
//  3. for each cell center, cut plane for is_insider(STL) (we need to play...)
//   3b. some cell centers rejected based on density per avg.neighbor-normal-distance-mm/min-voro-cell-mm

#include "common/vtk_common.h"
#include <voro++.hh>
#include <vtkSelectEnclosedPoints.h>
#include <vtkPolyData.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkCellLocator.h>
#include <vtkCellData.h>
#include <vtkSTLReader.h>
#include <cstdio>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>

using namespace voro;



//float MIN_CHANNEL_MM = -1.0;     // can't reduce it any further past this channel (+ 2x [extrude in mm])
float TUBE_RADIUS = -1.0;

const double LENGTH = 1.0;
const double WIDTH = 1.16; // 1.16" gives full diameter 1.75"
// length & width
const double height_multiplier = 2.54;
double total_length0 = LENGTH * height_multiplier;
const double width_multiplier = 2.54/2.0;
double total_width0 = WIDTH * width_multiplier;
// CYL voronoi zoom:
const double voronoi_zoom = 1.0;
double total_length = total_length0 * voronoi_zoom;
double total_width = total_width0 * voronoi_zoom;
// Set up constants for the container geometry
double xy_resolution = total_width;    // sets max x,y
const double x_min=-xy_resolution,x_max=xy_resolution;
const double y_min=-xy_resolution,y_max=xy_resolution;

double length_base = total_width;

const double z_min=length_base-total_length,z_max=length_base;
// Set the computational grid size
const int n_x=17,n_y=17,n_z=25;

float Gcell_face_mm;
// globals...
vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
vtkSmartPointer<vtkPolyData> stl_inner_shell = vtkSmartPointer<vtkPolyData>::New();
model_info info;	// TODO: maybe extrude_polydata_along_normal() can calc this itself, so don't need to pass to it...

/////////////////////////////////
////
//// Voro++ wall callbacks
////
/////////////////////////////////

/** \brief A class representing a user defiend wall.
 *
 * This class represents a shape defined by user callbacks. */

class wall_user_callback : public wall {
public:
	// return True if the point is inside, false if the point is outside.
	typedef bool(*myPointInsideCallback)(double, double, double);
	// return True if the cell still exists, false if the cell is deleted.
	typedef struct {
		double max_radius_squared;
		int number_of_faces;
		int number_of_edges;
		double total_edge_distance;
		double surface_area;
		double centroidx, centroidy, centroidz;
	} cell_info;
	typedef bool(*myCutCellCallback)(double, double, double, double*, double*, double*, double*, cell_info* ci);

	wall_user_callback(double xc_, double yc_, double zc_, myCutCellCallback cb_cp_, myPointInsideCallback cb_pi_ = NULL, int w_id_ = -99)
		: w_id(w_id_), xcb_point_inside(cb_pi_), xcb_cut_cell(cb_cp_), xc(xc_), yc(yc_), zc(zc_) {}
	bool point_inside(double x, double y, double z) {
		double t = xc + yc + zc;   // TODO: quiet compiler "warning: private field 'xc' is not used [-Wunused-private-field]"
		if (xcb_point_inside) return xcb_point_inside(x, y, z);
		return false;
	}
	/** Cuts a cell by the sphere wall object. The spherical wall is approximated by
	 * a single plane applied at the point on the sphere which is closest to the center
	 * of the cell. This works well for particle arrangements that are packed against
	 * the wall, but loses accuracy for sparse particle distributions.
	 * \param[in,out] c the Voronoi cell to be cut.
	 * \param[in] (x,y,z) the location of the Voronoi cell.
	 * \return True if the cell still exists, false if the cell is deleted. */
	template<class v_cell>
	bool cut_cell_base(v_cell &c, double x, double y, double z) {
		double xd, yd, zd, dq; xd = yd = zd = dq = 0.0;
		cell_info ci;
		ci.max_radius_squared = c.max_radius_squared();
		ci.number_of_faces = c.number_of_faces();
		ci.number_of_edges = c.number_of_edges();
		ci.total_edge_distance = c.total_edge_distance();
		ci.surface_area = c.surface_area();
		double d, f, g;
		c.centroid(d, f, g);
		ci.centroidx = d;
		ci.centroidy = f;
		ci.centroidz = g;
		if (xcb_cut_cell && xcb_cut_cell(x, y, z, &xd, &yd, &zd, &dq, &ci)) {
			//    if (dq>1e-5) {
			//        dq=2*(sqrt(dq)*rc-dq);
			return c.nplane(xd, yd, zd, dq, w_id);
		}
		return false;
	}

	bool cut_cell(voronoicell &c, double x, double y, double z) { return cut_cell_base(c, x, y, z); }
	bool cut_cell(voronoicell_neighbor &c, double x, double y, double z) { return cut_cell_base(c, x, y, z); }
private:
	const int w_id;
	myPointInsideCallback xcb_point_inside;
	myCutCellCallback xcb_cut_cell;
	const double xc, yc, zc;
};


// TODO: return false if not implemented yet
bool GetCellNormals(vtkPolyData* polydata, vtkIdType cellId, double* cN)
{
	//    std::cout << "Looking for cell normals..." << std::endl;
	//
	//    // Count points
	//    vtkIdType numCells = polydata->GetNumberOfCells();
	//    std::cout << "There are " << numCells << " cells." << std::endl;
	//
	//    // Count triangles
	//    vtkIdType numPolys = polydata->GetNumberOfPolys();
	//    std::cout << "There are " << numPolys << " polys." << std::endl;

	////////////////////////////////////////////////////////////////
	// Double normals in an array
	vtkDoubleArray* normalDataDouble = vtkDoubleArray::SafeDownCast(polydata->GetCellData()->GetArray("Normals"));

	if (normalDataDouble)
	{
		int nc = normalDataDouble->GetNumberOfTuples();
		std::cout << "There are " << nc
			<< " components in normalDataDouble" << std::endl;
		return false;//true;
	}

	////////////////////////////////////////////////////////////////
	// Double normals in an array
	vtkFloatArray* normalDataFloat = vtkFloatArray::SafeDownCast(polydata->GetCellData()->GetArray("Normals"));

	if (normalDataFloat)
	{
		int nc = normalDataFloat->GetNumberOfTuples();
		//        std::cout << "There are " << nc << " components in normalDataFloat" << std::endl;
		normalDataFloat->GetTuple(cellId, cN);
		//        cout << "Cell normal " << cellId << ": " << cN[0] << " " << cN[1] << " " << cN[2] << endl;
		return true;
	}

	////////////////////////////////////////////////////////////////
	// Point normals
	vtkDoubleArray* normalsDouble = vtkDoubleArray::SafeDownCast(polydata->GetCellData()->GetNormals());

	if (normalsDouble)
	{
		std::cout << "There are " << normalsDouble->GetNumberOfComponents()
			<< " components in normalsDouble" << std::endl;
		return false;//true;
	}

	////////////////////////////////////////////////////////////////
	// Point normals
	vtkFloatArray* normalsFloat = vtkFloatArray::SafeDownCast(polydata->GetCellData()->GetNormals());

	if (normalsFloat)
	{
		std::cout << "There are " << normalsFloat->GetNumberOfComponents()
			<< " components in normalsFloat" << std::endl;
		return false;//true;
	}

	/////////////////////////////////////////////////////////////////////
	// Generic type point normals
	vtkDataArray* normalsGeneric = polydata->GetCellData()->GetNormals(); //works
	if (normalsGeneric)
	{
		std::cout << "There are " << normalsGeneric->GetNumberOfTuples()
			<< " normals in normalsGeneric" << std::endl;

		double testDouble[3];
		normalsGeneric->GetTuple(0, testDouble);

		std::cout << "Double: " << testDouble[0] << " "
			<< testDouble[1] << " " << testDouble[2] << std::endl;

		// Can't do this:
		/*
		 float testFloat[3];
		 normalsGeneric->GetTuple(0, testFloat);

		 std::cout << "Float: " << testFloat[0] << " "
		 << testFloat[1] << " " << testFloat[2] << std::endl;
		 */
		return false;//true;
	}


	// If the function has not yet quit, there were none of these types of normals
	std::cout << "Normals not found!" << std::endl;
	return false;

}

// return True if the point is inside STL model, false if the point is outside.
static bool voro_cb_point_inside(double x, double y, double z) {

	double pt[3] = { x,y,z };
	bool is_inside = selectEnclosedPoints->IsInsideSurface(pt);
	bool is_inside_ex = false;
	printf("# VORO: voro_cb_point_inside:%0.8f %0.8f %0.8f:%0.8f %0.8f %0.8f: inside? ex%d:%d\n", pt[0], pt[1], pt[2], x, y, z, is_inside_ex, is_inside);

	// 031117: testing with slightly larger outside was not what wanted, voro needs to fit in given shape, so just look at those points
	return is_inside;// || is_inside_ex;
}

// return True if the cell still exists, false if the cell is deleted.
static bool voro_cb_cut_cell(double x, double y, double z,
	double* dx, double* dy, double* dz, double* dq, wall_user_callback::cell_info* ci) {

	// TODO: filter out some points when neighbor normals distance is increased
	if (!voro_cb_point_inside(x, y, z)) return false;

	double pt[3] = { x,y,z };
	//Find the closest points to TestPoint 'pt'
	double closestPoint[3];     //the coordinates of the closest point will be returned here
	double closestPointDist2;   //the squared distance to the closest point will be returned here
	vtkIdType cellId;           //the cell id of the cell containing the closest point will be returned here
	int subId;                  //this is rarely used (in triangle strips only, I believe)
	cellLocator->FindClosestPoint(pt, closestPoint, cellId, subId, closestPointDist2);
	//    printf("# VORO: voro_cb_cut_cell(%f,%f,%f) centroid(%f,%f,%f)\n",pt[0],pt[1],pt[2],center[0],center[1],center[2]);
	//    std::cout << "# VTK Coordinates of closest point: " << closestPoint[0] << " " << closestPoint[1] << " " << closestPoint[2] << std::endl;
	//    std::cout << "# VTK Squared distance to closest point: " << closestPointDist2 << std::endl;
	//    std::cout << "# VTK CellId: " << cellId << " subId: " << subId << std::endl;

	double dq_ret = closestPointDist2;
	// if too close to edge, say it's outside surface, one cell face_mm seems to work well, if you use Gcell_face/2.0, parts start sticking off
	float Gcell_face_mm2 = Gcell_face_mm * Gcell_face_mm;
	if (closestPointDist2 <= Gcell_face_mm2) {
		// todo: pull (x,y,z) back into surface by the amount (Gcell_face_mm - closestPointDist2)

//        return false;
		dq_ret = Gcell_face_mm2 - (Gcell_face_mm2 - closestPointDist2);
		//printf("# VORO:voro_cb_cut_cell(%f,%f,%f) is close to surface (%f mm^2), set *dq=%f\n",x,y,z,closestPointDist2,dq_ret);

	}

	// Try to read normals
	double cN[3];
	if (!GetCellNormals(stl_inner_shell, cellId, cN))
	{
		std::cout << "ERROR: No cell normals were found" << std::endl;
		return false;
	}
	//    printf("# VTK NORMAL:%f,%f,%f\n",cN[0],cN[1],cN[2]);
	//    printf("# VORO:radius^2:%f faces:%d edges:%d edge_dist:%f surface_area:%f\n",
	//           ci->max_radius_squared,ci->number_of_faces,ci->number_of_edges,ci->total_edge_distance,ci->surface_area);
	//    //    printf("# main:: voro_cb_cut_cell:%f,%f,%f: inside? 1\n",pt[0],pt[1],pt[2]);
	*dx = cN[0];
	*dy = cN[1];
	*dz = cN[2];
	*dq = dq_ret;
	//    printf("# VORODIST:%f = %f\n",closestPointDist2,dq_ret);
	return true;
}

void config_enclosedPoints_inner(double scale) {

#if 0
	//    vtkSmartPointer<vtkPolyData> temppolydata = vtkSmartPointer<vtkPolyData>::New();
	//    temppolydata = reader->GetOutput();
	//
	//    // Create the locator
	//    vtkSmartPointer<vtkOBBTree> tree = vtkSmartPointer<vtkOBBTree>::New();
	//    tree->SetDataSet(temppolydata);
	//    tree->BuildLocator();
	//    
	//    // call fcn, voro_fill also needs this
	//    int ret = -1;
	//    vtkSmartPointer<vtkPolyData> temppolydata2 = vtkSmartPointer<vtkPolyData>::New();
	//    temppolydata2 = extrude_polydata_along_normal(temppolydata,scale,tree,MIN_CHANNEL_MM,true,info,&ret);
	//    if (ret != 0) {
	//        printf("ERROR: error returned from extrude_polydata_along_normal()\n");
	//        return;// ret;
	//    }
#else
	vtkSmartPointer<vtkPolyData> temppolydata2 = vtkSmartPointer<vtkPolyData>::New();
	temppolydata2 = reader->GetOutput();
#endif

	// for the inner surface
	vtkSmartPointer<vtkPolyDataNormals> pnormal = vtkSmartPointer<vtkPolyDataNormals>::New();
#if VTK_MAJOR_VERSION <= 5
	pnormal->SetInput(temppolydata2);
#else
	pnormal->SetInputData(temppolydata2);
#endif
	pnormal->SplittingOff();
	pnormal->FlipNormalsOff();
	pnormal->ComputeCellNormalsOn();
	pnormal->ComputePointNormalsOff();
	pnormal->Update();

	stl_inner_shell = pnormal->GetOutput();
	selectEnclosedPoints->Initialize(stl_inner_shell);
	selectEnclosedPoints->SetTolerance(0.0000001);
	//    selectEnclosedPoints->Update();
}


int main(int argc, char *argv[])
{

	// open STL
	//std::string inputFilename = "UMS5_bunny.stl";
	std::string inputFilename = "UMS5_Knee_Prosthesis.stl";
	read_stl(&reader, inputFilename, true);

	//Gcell_face_mm = 1;

	get_polydata_info(&info, reader->GetOutput(), true);
	
	// config enclosed points
	config_enclosedPoints_inner(1);

	// INIT: closest points in reader() STL data
	cellLocator->SetDataSet(stl_inner_shell);
	cellLocator->BuildLocator();
	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block.
    container con(info.min.x,info.max.x,info.min.y,info.max.y,info.min.z,info.max.z,n_x,n_y,n_z,
                  false,false,false,8);
    
	
    wall_user_callback wall(info.center[0],info.center[1],info.center[2],&voro_cb_cut_cell,&voro_cb_point_inside);
    con.add_wall(wall);
    con.import("pack_2");
    con.draw_cells_pov("test_v.pov");
	//con.print_custom("%B","cad_small-cylinder.csv");
	con.draw_particles_pov("test_p.pov");
    
    return 0;
}











