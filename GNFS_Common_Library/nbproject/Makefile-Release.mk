#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=
AS=

# Macros
CND_PLATFORM=GNU-MacOSX
CND_CONF=Release
CND_DISTDIR=dist

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=build/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/Sieving.o \
	${OBJECTDIR}/LineSieving.o \
	${OBJECTDIR}/gnfs_utils.o \
	${OBJECTDIR}/Relation.o \
	${OBJECTDIR}/pair_intType_intType.o \
	${OBJECTDIR}/triple_intType.o \
	${OBJECTDIR}/FactorBases.o \
	${OBJECTDIR}/LinearAlgebra.o \
	${OBJECTDIR}/SquareRoot.o \
	${OBJECTDIR}/pair_long_long.o \
	${OBJECTDIR}/Parameters.o \
	${OBJECTDIR}/PolynomialSelection.o \
	${OBJECTDIR}/timer.o \
	${OBJECTDIR}/vec_intType.o \
	${OBJECTDIR}/MatrixConstruction.o \
	${OBJECTDIR}/Initialization.o \
	${OBJECTDIR}/pair_uintType_uintType.o \
	${OBJECTDIR}/vec_uintType.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-m64
CXXFLAGS=-m64

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	${MAKE}  -f nbproject/Makefile-Release.mk dist/Release/libgnfs_common_library.a

dist/Release/libgnfs_common_library.a: ${OBJECTFILES}
	${MKDIR} -p dist/Release
	${RM} dist/Release/libgnfs_common_library.a
	${AR} rv dist/Release/libgnfs_common_library.a ${OBJECTFILES} 
	$(RANLIB) dist/Release/libgnfs_common_library.a

${OBJECTDIR}/Sieving.o: nbproject/Makefile-${CND_CONF}.mk Sieving.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/Sieving.o Sieving.cpp

${OBJECTDIR}/LineSieving.o: nbproject/Makefile-${CND_CONF}.mk LineSieving.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/LineSieving.o LineSieving.cpp

${OBJECTDIR}/gnfs_utils.o: nbproject/Makefile-${CND_CONF}.mk gnfs_utils.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/gnfs_utils.o gnfs_utils.cpp

${OBJECTDIR}/Relation.o: nbproject/Makefile-${CND_CONF}.mk Relation.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/Relation.o Relation.cpp

${OBJECTDIR}/pair_intType_intType.o: nbproject/Makefile-${CND_CONF}.mk pair_intType_intType.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/pair_intType_intType.o pair_intType_intType.cpp

${OBJECTDIR}/triple_intType.o: nbproject/Makefile-${CND_CONF}.mk triple_intType.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/triple_intType.o triple_intType.cpp

${OBJECTDIR}/FactorBases.o: nbproject/Makefile-${CND_CONF}.mk FactorBases.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/FactorBases.o FactorBases.cpp

${OBJECTDIR}/LinearAlgebra.o: nbproject/Makefile-${CND_CONF}.mk LinearAlgebra.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/LinearAlgebra.o LinearAlgebra.cpp

${OBJECTDIR}/SquareRoot.o: nbproject/Makefile-${CND_CONF}.mk SquareRoot.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/SquareRoot.o SquareRoot.cpp

${OBJECTDIR}/pair_long_long.o: nbproject/Makefile-${CND_CONF}.mk pair_long_long.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/pair_long_long.o pair_long_long.cpp

${OBJECTDIR}/Parameters.o: nbproject/Makefile-${CND_CONF}.mk Parameters.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/Parameters.o Parameters.cpp

${OBJECTDIR}/PolynomialSelection.o: nbproject/Makefile-${CND_CONF}.mk PolynomialSelection.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/PolynomialSelection.o PolynomialSelection.cpp

${OBJECTDIR}/timer.o: nbproject/Makefile-${CND_CONF}.mk timer.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/timer.o timer.cpp

${OBJECTDIR}/vec_intType.o: nbproject/Makefile-${CND_CONF}.mk vec_intType.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/vec_intType.o vec_intType.cpp

${OBJECTDIR}/MatrixConstruction.o: nbproject/Makefile-${CND_CONF}.mk MatrixConstruction.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/MatrixConstruction.o MatrixConstruction.cpp

${OBJECTDIR}/Initialization.o: nbproject/Makefile-${CND_CONF}.mk Initialization.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/Initialization.o Initialization.cpp

${OBJECTDIR}/pair_uintType_uintType.o: nbproject/Makefile-${CND_CONF}.mk pair_uintType_uintType.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/pair_uintType_uintType.o pair_uintType_uintType.cpp

${OBJECTDIR}/vec_uintType.o: nbproject/Makefile-${CND_CONF}.mk vec_uintType.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -O2 -Iinclude -MMD -MP -MF $@.d -o ${OBJECTDIR}/vec_uintType.o vec_uintType.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	${RM} -r build/Release
	${RM} dist/Release/libgnfs_common_library.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
