# =========================================================================
# Compiler and Flags Configuration
# =========================================================================
#FC       =   nvfortran
#FCFLAGS  = -O3 -tp=native #-ffree-line-length-none
FC       =   gfortran -ffast-math -funroll-loops -flto # -fcheck-array-temporaries  -fopt-info-vec-missed 
FCFLAGS  = -O3 -march=skylake-avx512 -ffree-line-length-none  #  -g    -fbounds-check

#FC       =   ifort
#FCFLAGS  =    -O3      -xCORE-AVX512 -qopt-zmm-usage=high  #-traceback -O2 # -CB  -g -traceback #-Wall -Wextra -fcheck=all

# gfortran option -fcheck=bounds -g -fbacktrace
# =========================================================================
# Files and Targets
# =========================================================================
TARGET   = MMM_1Dconvect

# Explicit source listing tracking compilation dependencies
SRCS     = module_SchemeSetup.f90\
           module_VTK_MultiMoment_Writer.f90\
           module_froto_vortex_functions.f90 \
           MMM_module_spatial.f90 \
           MMM_module_TimeINT.f90 \
           MMM_1Dconvect.f90

# Maps source files to their corresponding object file names
OBJS     = $(SRCS:.f90=.o)

# =========================================================================
# Build Rules
# =========================================================================

# Default target rule
all: $(TARGET)

# Linking rule to create the final executable binary
$(TARGET): $(OBJS)
	$(FC) $(FCFLAGS) -o $(TARGET) $(OBJS)

# General pattern rule for compiling Fortran source files into object files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -o $@

# =========================================================================
# Explicit Fortran Module Dependencies
# =========================================================================
# 1. module_froto_vortex_functions is entirely independent
module_VTK_MultiMoment_Writer.o: module_VTK_MultiMoment_Writer.f90

module_SchemeSetup.o: module_SchemeSetup.f90

module_froto_vortex_functions.o: module_froto_vortex_functions.f90 module_SchemeSetup.o

# 2. MMM_module_spatial relies on module_froto_vortex_functions.o
MMM_module_spatial.o: MMM_module_spatial.f90 module_froto_vortex_functions.o module_SchemeSetup.o

# 3. MMM_module_TimeINT uses MMM_module_spatial.mod
MMM_module_TimeINT.o: MMM_module_TimeINT.f90 MMM_module_spatial.o  

# 4. MMM_1Dconvect uses MMM_module_TimeINT.mod
MMM_1Dconvect.o: MMM_1Dconvect.f90 MMM_module_TimeINT.o module_VTK_MultiMoment_Writer.o MMM_module_TimeINT.o

# =========================================================================
# Maintenance Operations
# =========================================================================
.PHONY: clean run

# Cleans up binary files, object files, and compiled Fortran module interfaces
clean:
	rm -f $(TARGET) $(OBJS) *.mod  

# Convenience target to automatically recompile and run the advection simulation
run: $(TARGET)
	./$(TARGET)
	
	
