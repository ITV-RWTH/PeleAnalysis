- __Step 1:__ ```cd $(PELE_ANALYSIS_HOME)/Submodules/PelePhysics/Submodules/sundials```

- __Step 2:__ ```mkdir build && cd build```

- __Step 3:__ ```cmake .. \
                  -DCMAKE_INSTALL_PREFIX=/pfad/zu/sundials-install \
                  -DENABLE_MPI=ON \
                  -DSUNDIALS_INDEX_SIZE=32   #IMPORTANT!```

- __Step 4:__ ```make -j && make install```

- __Step 5:__ Make sure that ```SUNDIALS_HOME``` in ```$(PELE_ANALYSIS_HOME)/Tools/GNUmake/Make.ModelSpecific``` points to ```$(PELE_PHYSICS_HOME)/Submodules/sundials/build/install``` (the ```/build/install```-part is absolutely crucial and may be missing!)
