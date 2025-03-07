#!/bin/sh
skip=49

tab='	'
nl='
'
IFS=" $tab$nl"

umask=`umask`
umask 77

gztmpdir=
trap 'res=$?
  test -n "$gztmpdir" && rm -fr "$gztmpdir"
  (exit $res); exit $res
' 0 1 2 3 5 10 13 15

case $TMPDIR in
  / | /*/) ;;
  /*) TMPDIR=$TMPDIR/;;
  *) TMPDIR=/tmp/;;
esac
if type mktemp >/dev/null 2>&1; then
  gztmpdir=`mktemp -d "${TMPDIR}gztmpXXXXXXXXX"`
else
  gztmpdir=${TMPDIR}gztmp$$; mkdir $gztmpdir
fi || { (exit 127); exit 127; }

gztmp=$gztmpdir/$0
case $0 in
-* | */*'
') mkdir -p "$gztmp" && rm -r "$gztmp";;
*/*) gztmp=$gztmpdir/`basename "$0"`;;
esac || { (exit 127); exit 127; }

case `printf 'X\n' | tail -n +1 2>/dev/null` in
X) tail_n=-n;;
*) tail_n=;;
esac
if tail $tail_n +$skip <"$0" | gzip -cd > "$gztmp"; then
  umask $umask
  chmod 700 "$gztmp"
  (sleep 5; rm -fr "$gztmpdir") 2>/dev/null &
  "$gztmp" ${1+"$@"}; res=$?
else
  printf >&2 '%s\n' "Cannot decompress $0"
  (exit 127); res=127
fi; exit $res
�,j�gsetup.sh �Z{sU�{�S�֞�W�5�PK�<
�ܭ����;3������A��"���V-�SEE��$�+� ~���L��+�}�9/��hRd��=��s�9��߹ö��YUg%��q��{���{�s%��u�^�\��Gz��O���Z\��eG5t�wUx���LT�v$Mٵ4�Yr�)"���:��P��?��Txi��Jx݆��vI-i.SU��;��Ē�!������ϒ��*u�ϮHY��v�%�3vm4`>W%AtC�kD�2����z�_&�T�^%�����s�$D	%:���xx�y�kC����C�`�Wf��lz�J����I�YE���.O��ʗ�b�6H�B��M���"��	��MÅ9Iw�1"�PJT�8�^�Mb�MQ7D���W/w�o��h	x�6]�[�������6>���
LM�ݟ�\9���]�<LO��Ftݠ��u�#:����݀k5I�ʺh�N�B H�kD�BS���u�p]�5�����lvM�L�Kvb�>�Ӎ^��0������t۠��B��o֮-�o�ſ����w�Onyw.{�o��T��뭕�l�[��-]�VF<��ρ� /Gm��F��vuE����#��=v9�0�C�q�j�f߷�\Ô�Gg߾	�"�lgŘ�5CR��LV�OI�mx|q���lf�!(IF�c���}r����?��$J���ޙ+��7�O�]=�Y��xK����{����lk�G.e@?�2�kx�G�ॗ`��>n
��0�i��IQ,(��W����
�\f�y"�6�2���1�D�E����w�.�בֿ|��qS��mצ9�i����i�޵���\��Ls�`�E��si	�rc�<G��c]�1n�-���ڗ˘��i���b^��ؾ�d��Cf��������M�*��~�����=�X�2X�<x�Tk僟]MR����l҇eL��������©�Ə�|��%g��U���&.i(�'��`���Y����LsĖ-դ���u�������e�8s�U�`~���֌j�ƽ����"��8Ls�0�˶�05�M��ϒ�-R�ᏭbB!O����*�}���~�c�}�D��,7\�QZlB��ofv��@���&+�P���.�j�^Djۯ�EJǵS���� e��y���K7���Pt��g#0�%~�K�n����̫� p�dJ֝�����U"�B��^K� h�.T-��W�Gx�er\�]M��l�@� #9]�h�[�Ҿ�J�H��#�f3�Hv�y���nF�Q�NčK�v�r����4MY�� Ȉyk�N䅄#����gvi=�h8tk�J�E��$�1�#��W?��1���*��۹����0넜�@us8�i *C'�@t�*C��z�I$?����Y�K�{���g��T�<=\������� �D��:w�Z���v�v{�>�'&��z��y�߳��]�=8�������_B�����rl�����s�ڊ��=qQU����k���:�'�4����c�G����'�lX��lVC]j���YS�6�-��|>�ϿY߉3X��f4?WSK�	���Q��۬�Q{��Z�4��]�?H� m(�%��U��w� ���}�|)����Ӟ��EN��F40�YM���R�g�!/��^g%�{MI�(O5�QFnlwi3�]�<(�2�GOf�$$������%�Ϸ�!�sև�Y��?�� g�y�'75Tojb��o%���m�ʬxgզ�G0�t�ȋ�b�>�4�-_m=8��ft;��zONs ס�^UHEB�mn3�z�x:�ӿ)��GЦ�s`�����Fy:� IIf�a�:b3 -rQ_~H"�h��'E>�=�ei���t��A"���y��U��y7��0D��� [���-��b	b�I������B��Y�e4�W��̹�!�	=Ḥ�$��&�Z�
9\2�.���V�;ŇӴ���L�����,	��'�G<�S��������:�U�7X��G</��	$�A�\�Dw�����E�I�GHU��� ��H6\��¿��R1dA�#JX�|��=̱��o�.Ao%.D3W}bz_�W�ꊤ:��_�-�G}�2�J���Vp@�^�WT��H=f��Y��2H�X]��.�pbS<d���������Mr=D�$���ݟ*>AF5b.l�a�g01l5����I���IM;�?0g��X,��},_�]�DW1�Ǉ۠$J�)R$���'��ym�ȱ��&&�Iw%��`|b�Hٗ�=3q������>�U�Z����&����J����*�M�4����Ж�Yh�����~��]Ayـ2��}�3z��l�T}��������~��;�[����{�ޙ�ąɮs3DcVsƵ%lJ�5���]ݥW,�ɫ��3���յ���U������=�t7)C�+[(L����?�q���whc���w�����[�Χ8��Vu��a�����ʷ�W�Mn4����e��c�$:��ZI�|L���%�^&I
1I�r���x:�<��ب�+82~ ��HNfG�F�g.Gg{O�]��.������k�<���v��6���h�p,���C�`vi@���S�1��zؑ������,���P����h�k�yWM�����Ӕh�A��4et�:���&O��>�ci➳9h]v���U��ه��W�׿K���}@3���ُ��#G��v�	�|�~��o_��ݹ�#1�OlI����������/,w~��]�<
�L&�'BB��_�����M.#  