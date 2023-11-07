set(TEST_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})

configure_file(caradonna_finite_span_subsonic_entropy_grid_adaptation.prm caradonna_finite_span_subsonic_entropy_grid_adaptation.prm COPYONLY)
add_test(
	NAME 3D_EULER_CARADONNA_FINITE_SUBSONIC_ENTROPY_ADAPTATION
	COMMAND mpirun -np ${MPIMAX} ${EXECUTABLE_OUTPUT_PATH}/PHiLiP_3D -i ${CMAKE_CURRENT_BINARY_DIR}/caradonna_finite_span_subsonic_entropy_grid_adaptation.prm
	WORKING_DIRECTORY ${TEST_OUTPUT_DIR}
)

configure_file(caradonna_finite_span_transonic_entropy_grid_adaptation.prm caradonna_finite_span_transonic_entropy_grid_adaptation.prm COPYONLY)
add_test(
	NAME 3D_EULER_CARADONNA_FINITE_TRANSONIC_ENTROPY_ADAPTATION
	COMMAND mpirun -np ${MPIMAX} ${EXECUTABLE_OUTPUT_PATH}/PHiLiP_3D -i ${CMAKE_CURRENT_BINARY_DIR}/caradonna_finite_span_transonic_entropy_grid_adaptation.prm
	WORKING_DIRECTORY ${TEST_OUTPUT_DIR}
)

configure_file(caradonna_infinite_span_subsonic_entropy_grid_adaptation.prm caradonna_infinite_span_subsonic_entropy_grid_adaptation.prm COPYONLY)
add_test(
	NAME 3D_EULER_CARADONNA_INFINITE_SUBSONIC_ENTROPY_ADAPTATION
	COMMAND mpirun -np ${MPIMAX} ${EXECUTABLE_OUTPUT_PATH}/PHiLiP_3D -i ${CMAKE_CURRENT_BINARY_DIR}/caradonna_infinite_span_subsonic_entropy_grid_adaptation.prm
	WORKING_DIRECTORY ${TEST_OUTPUT_DIR}
)

configure_file(caradonna_infinite_span_transonic_entropy_grid_adaptation.prm caradonna_infinite_span_transonic_entropy_grid_adaptation.prm COPYONLY)
add_test(
	NAME 3D_EULER_CARADONNA_INFINITE_TRANSONIC_ENTROPY_ADAPTATION
	COMMAND mpirun -np ${MPIMAX} ${EXECUTABLE_OUTPUT_PATH}/PHiLiP_3D -i ${CMAKE_CURRENT_BINARY_DIR}/caradonna_infinite_span_transonic_entropy_grid_adaptation.prm
	WORKING_DIRECTORY ${TEST_OUTPUT_DIR}
)