sed \
	-e 's/#define *MS_CL_COMPILER/#undef    MS_CL_COMPILER/' \
	-e 's/#define *GNU_GPP_COMPILER/#undef    GNU_GPP_COMPILER/' \
	-e 's/#define *SUN_CC_COMPILER/#undef    SUN_CC_COMPILER/' \
	-e 's/#define *PGI_PGCC_COMPILER/#undef    PGI_PGCC_COMPILER/' \
	scltypes_float.tpl > scltypes.tmp

sed -e 's/#undef *MS_CL_COMPILER/#define   MS_CL_COMPILER/' \
	scltypes.tmp > scltypes.ms
sed -e 's/#undef *GNU_GPP_COMPILER/#define   GNU_GPP_COMPILER/' \
	scltypes.tmp > scltypes.gpp
sed -e 's/#undef *SUN_CC_COMPILER/#define   SUN_CC_COMPILER/' \
	scltypes.tmp > scltypes.sun
sed -e 's/#undef *PGI_PGCC_COMPILER/#define   PGI_PGCC_COMPILER/' \
	scltypes.tmp > scltypes.pgi


echo " "
echo "scltypes.ms:"
grep COMPILER scltypes.ms

echo " "
echo "scltypes.gpp:"
grep COMPILER scltypes.gpp

echo " "
echo "scltypes.sun:"
grep COMPILER scltypes.sun

echo " "
echo "scltypes.pgi:"
grep COMPILER scltypes.pgi

rm -f scltypes.tmp

echo " "
