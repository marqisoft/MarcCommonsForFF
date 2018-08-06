### Instructions for setting up a Visual Studio FEFLOW-IFM project to locate MarcCommonsForFF and netCDF library source codes
*(e.g., the PaleoSea2D plugin's project)*

1. Open the C++ solution (project) in Microsoft Visual Studio.

2. Open the project's properties from the menu (Project > ... Properties).

3. In the "C/C++" group > "General" subgroup:
   click on the text field of Additional Include directories, then\
   click on the down-arrow button (to the right), then\
   click "Edit...", and finally:\
   add the two paths below*, after the default 'ifm' paths already present:\
   C:\Users\USERNAME\Documents\Visual Studio 2015\Projects\MarcCommonsForFF\
   C:\Program Files\netCDF441\include\
   (\* Adapt, of course, the folder paths to your own custom setup.)

4. In the "Linker" group > "General" subgroup:
   click on the text field of Additional Library Directories, then\
   click on the down-arrow button (to the right), then\
   click "Edit...", and finally:\
   add the path below* (normally will be the only one in the multiline text zone)\
   C:\Program Files\netCDF441\lib\
   (\* Adapt, of course and again, the folder paths to your own custom setup.)

5. In the "Linker" group > "Input" subgroup:
   click on the text field of Additional Dependencies, then\
   click on the down-arrow button (to the right), then\
   click "Edit...", and finally:\
   add the filename below (normally will be the only one in the multiline text zone)\
   netcdf.lib

These settings should be sufficient for the interpreter & compiler to locate and use
the source code files and libraries required for instance by the PaleoSea2D plugin.

Copyright (c) 2018 Marc Laurencelle (alias marQIsoft)
