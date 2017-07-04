CFLAGS = -I${PETSC_DIR}/include/

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

HeatEqnFDM_Petsc: HeatEqnFDM_Petsc.o chkopts
	-${CLINKER} -o HeatEqnFDM_Petsc HeatEqnFDM_Petsc.o  ${PETSC_SYS_LIB}
	${RM} -f HeatEqnFDM_Petsc.o

