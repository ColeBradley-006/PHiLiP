set(TEST_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})

configure_file(caradonna_finite.prm caradonna_finite.prm COPYONLY)
add_test(
	NAME 3D_EULER_CARADONNA_FINITE
	COMMAND mpirun -np ${MPIMAX} ${EXECUTABLE_OUTPUT_PATH}/PHiLiP_3D -i ${CMAKE_CURRENT_BINARY_DIR}/caradonna_finite.prm
	WORKING_DIRECTORY ${TEST_OUTPUT_DIR}
)

configure_file(caradonna_naca.prm caradonna_naca.prm COPYONLY)
add_test(
	NAME 3D_EULER_CARADONNA_NACA
	COMMAND mpirun -np ${MPIMAX} ${EXECUTABLE_OUTPUT_PATH}/PHiLiP_3D -i ${CMAKE_CURRENT_BINARY_DIR}/caradonna_naca.prm
	WORKING_DIRECTORY ${TEST_OUTPUT_DIR}
)

configure_file(caradonna_span.prm caradonna_span.prm COPYONLY)
add_test(
	NAME 3D_EULER_CARADONNA_SPAN
	COMMAND mpirun -np ${MPIMAX} ${EXECUTABLE_OUTPUT_PATH}/PHiLiP_3D -i ${CMAKE_CURRENT_BINARY_DIR}/caradonna_span.prm
	WORKING_DIRECTORY ${TEST_OUTPUT_DIR}
)