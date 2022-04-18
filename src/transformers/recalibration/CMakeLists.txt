

file(GLOB RECALIBRATION_SOURCE "./*.cpp")
add_library(Recalibration SHARED ${RECALIBRATION_SOURCE}  ${PROJECT_SOURCE_DIR}/src/MathUtils.cpp)