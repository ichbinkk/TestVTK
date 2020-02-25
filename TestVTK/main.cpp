//#include <vtkAutoInit.h>
//#define vtkRenderingCore_AUTOINIT 2(vtkRenderingOpenGL2, vtkInteractionStyle)
//
//#include <vtkSphereSource.h>
//#include <vtkPolyData.h>
//#include <vtkSmartPointer.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindowInteractor.h>
//
//#include <vtkLineSource.h>
//#include <vtkTubeFilter.h>
//
//#include <vtkBooleanOperationPolyDataFilter.h>
//std::string operation("intersection");


//int main(int, char *[])
//{
//	// Create a sphere
//	vtkSmartPointer<vtkSphereSource> sphereSource =
//		vtkSmartPointer<vtkSphereSource>::New();
//	sphereSource->SetCenter(0.0, 0.0, 0.0);
//	sphereSource->SetRadius(5.0);
//
//	//mapper
//	vtkSmartPointer<vtkPolyDataMapper> mapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	mapper->SetInputConnection(sphereSource->GetOutputPort());
//
//	//actor
//	vtkSmartPointer<vtkActor> actor =
//		vtkSmartPointer<vtkActor>::New();
//	actor->SetMapper(mapper);
//
//	//renderer ,renderWindow, renderWindowInteractor.
//	vtkSmartPointer<vtkRenderer> renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	renderer->AddActor(actor);
//	renderer->SetBackground(.3, .6, .3); // Background color green
//	renderWindow->Render();
//	renderWindowInteractor->Start();
//	return EXIT_SUCCESS;
//}

//测试vtk的布尔操作
//int main(int argc, char *argv[])
//{
//	vtkSmartPointer<vtkPolyData> input1;
//	vtkSmartPointer<vtkPolyData> input2;
//
//	vtkSmartPointer<vtkSphereSource> sphereSource1 = vtkSmartPointer<vtkSphereSource>::New();
//	sphereSource1->SetCenter(.25, 0, 0);
//	sphereSource1->Update();
//	input1 = sphereSource1->GetOutput();// use the polydata as input
//
//	vtkSmartPointer<vtkSphereSource> sphereSource2 = vtkSmartPointer<vtkSphereSource>::New();
//	sphereSource2->Update();
//	input2 = sphereSource2->GetOutput();
//
//	vtkSmartPointer<vtkBooleanOperationPolyDataFilter> booleanOperation =
//		vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
//	if (operation == "union")
//	{
//		booleanOperation->SetOperationToUnion();
//	}
//	else if (operation == "intersection")
//	{
//		booleanOperation->SetOperationToIntersection();
//	}
//	else if (operation == "difference")
//	{
//		booleanOperation->SetOperationToDifference();
//	}
//	else
//	{
//		std::cout << "Unknown operation: " << operation << std::endl;
//		return EXIT_FAILURE;
//	}
//#if VTK_MAJOR_VERSION <= 5
//	booleanOperation->SetInputConnection(0, input1->GetProducerPort());
//	booleanOperation->SetInputConnection(1, input2->GetProducerPort());
//#else
//	booleanOperation->SetInputData(0, input1); // set the input data
//	booleanOperation->SetInputData(1, input2);
//#endif
//	vtkSmartPointer<vtkPolyDataMapper> booleanOperationMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	booleanOperationMapper->SetInputConnection(booleanOperation->GetOutputPort());
//	booleanOperationMapper->ScalarVisibilityOff();
//
//	vtkSmartPointer<vtkActor> booleanOperationActor =
//		vtkSmartPointer<vtkActor>::New();
//	booleanOperationActor->SetMapper(booleanOperationMapper);
//
//	vtkSmartPointer<vtkRenderer> renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	renderer->AddViewProp(booleanOperationActor);
//	renderer->SetBackground(.1, .2, .3);
//	vtkSmartPointer<vtkRenderWindow> renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//
//	vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renWinInteractor->SetRenderWindow(renderWindow);
//
//	renderWindow->Render();
//	renWinInteractor->Start();
//
//	return EXIT_SUCCESS;
//}

//测试点是否在vtkpolydata内部
//#include <vtkVersion.h>
//#include <vtkPolyData.h>
//#include <vtkPointData.h>
//#include <vtkCubeSource.h>
//#include <vtkSmartPointer.h>
//#include <vtkSelectEnclosedPoints.h>
//#include <vtkIntArray.h>
//#include <vtkDataArray.h>
//#include <vtkVertexGlyphFilter.h>
//#include <vtkProperty.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkActor.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkVertexGlyphFilter.h>
//
//
//int main(int, char *argv[])
//{
//	//cube centered in origin, 1cm side.
//	vtkSmartPointer<vtkCubeSource> cubeSource =
//		vtkSmartPointer<vtkCubeSource>::New();
//	cubeSource->Update();
//
//	vtkPolyData* cube = cubeSource->GetOutput();
//
//	double testInside[3] = { 0.0, 0.0, 0.0 };
//	double testOutside[3] = { 0.7, 0.0, 0.0 };
//	double testBorderOutside[3] = { 0.5, 0.0, 0.0 };
//	vtkSmartPointer<vtkPoints> points =
//		vtkSmartPointer<vtkPoints>::New();
//	points->InsertNextPoint(testInside);
//	points->InsertNextPoint(testOutside);
//	points->InsertNextPoint(testBorderOutside);
//
//	vtkSmartPointer<vtkPolyData> pointsPolydata =
//		vtkSmartPointer<vtkPolyData>::New();
//	pointsPolydata->SetPoints(points);
//
//	//Points inside test
//	vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints =
//		vtkSmartPointer<vtkSelectEnclosedPoints>::New();
//#if VTK_MAJOR_VERSION <= 5
//	selectEnclosedPoints->SetInput(pointsPolydata);
//#else
//	selectEnclosedPoints->SetInputData(pointsPolydata);
//#endif
//#if VTK_MAJOR_VERSION <= 5
//	selectEnclosedPoints->SetSurface(cube);
//#else
//	selectEnclosedPoints->SetSurfaceData(cube);
//#endif
//	selectEnclosedPoints->Update();
//
//	for (unsigned int i = 0; i < 3; i++)
//	{
//		std::cout << "Point " << i << ": " << selectEnclosedPoints->IsInside(i) << std::endl;
//	}
//
//	vtkDataArray* insideArray = vtkDataArray::SafeDownCast(selectEnclosedPoints->GetOutput()->GetPointData()->GetArray("SelectedPoints"));
//
//	for (vtkIdType i = 0; i < insideArray->GetNumberOfTuples(); i++)
//	{
//		std::cout << i << " : " << insideArray->GetComponent(i, 0) << std::endl;
//	}
//
//
//	//RENDERING PART
//
//	//Cube mapper, actor
//	vtkSmartPointer<vtkPolyDataMapper> cubeMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	cubeMapper->SetInputConnection(cubeSource->GetOutputPort());
//
//	vtkSmartPointer<vtkActor> cubeActor =
//		vtkSmartPointer<vtkActor>::New();
//	cubeActor->SetMapper(cubeMapper);
//	cubeActor->GetProperty()->SetOpacity(0.5);
//
//	//Points mapper, actor
//	//First, apply vtkVertexGlyphFilter to make cells around points, vtk only render cells.
//	vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter =
//		vtkSmartPointer<vtkVertexGlyphFilter>::New();
//#if VTK_MAJOR_VERSION <= 5
//	vertexGlyphFilter->AddInput(pointsPolydata);
//#else
//	vertexGlyphFilter->AddInputData(pointsPolydata);
//#endif
//	vertexGlyphFilter->Update();
//
//	vtkSmartPointer<vtkPolyDataMapper> pointsMapper =
//		vtkSmartPointer<vtkPolyDataMapper>::New();
//	pointsMapper->SetInputConnection(vertexGlyphFilter->GetOutputPort());
//
//	vtkSmartPointer<vtkActor> pointsActor =
//		vtkSmartPointer<vtkActor>::New();
//	pointsActor->SetMapper(pointsMapper);
//	pointsActor->GetProperty()->SetPointSize(5);//定义点的尺寸大小，这样点才能在画布上显示出来
//	pointsActor->GetProperty()->SetColor(0.0, 0.0, 1);
//
//	//Create a renderer, render window, and interactor
//	vtkSmartPointer<vtkRenderer> renderer =
//		vtkSmartPointer<vtkRenderer>::New();
//	vtkSmartPointer<vtkRenderWindow> renderWindow =
//		vtkSmartPointer<vtkRenderWindow>::New();
//	renderWindow->AddRenderer(renderer);
//	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//		vtkSmartPointer<vtkRenderWindowInteractor>::New();
//	renderWindowInteractor->SetRenderWindow(renderWindow);
//
//	// Add the actor to the scene
//	renderer->AddActor(cubeActor);
//	renderer->AddActor(pointsActor);
//	renderer->SetBackground(.0, 1, .0);
//
//	// Render and interact
//	renderWindow->SetWindowName(argv[0]);
//	renderWindow->Render();
//	renderWindowInteractor->Start();
//
//	return EXIT_SUCCESS;
//}







