gfortran ./Code/*.F90 ./lib/*.lib ^
         -o main.exe -J ./x64/gfortran/Release/ ^
         --debug -fbackslash -O3