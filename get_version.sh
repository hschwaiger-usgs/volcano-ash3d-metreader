echo -n "      character(len=40),parameter,public :: MR_GitComID ='" > MR_version.h
git log -n 1 | grep commit | cut -f 2 -d' ' | tr -d $'\n' >> MR_version.h
echo -n "'" >> MR_version.h
