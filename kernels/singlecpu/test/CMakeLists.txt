####################################################################
# Copyright © 2018 Computational Molecular Biology Group,          #
#                  Freie Universität Berlin (GER)                  #
#                                                                  #
# Redistribution and use in source and binary forms, with or       #
# without modification, are permitted provided that the            #
# following conditions are met:                                    #
#  1. Redistributions of source code must retain the above         #
#     copyright notice, this list of conditions and the            #
#     following disclaimer.                                        #
#  2. Redistributions in binary form must reproduce the above      #
#     copyright notice, this list of conditions and the following  #
#     disclaimer in the documentation and/or other materials       #
#     provided with the distribution.                              #
#  3. Neither the name of the copyright holder nor the names of    #
#     its contributors may be used to endorse or promote products  #
#     derived from this software without specific                  #
#     prior written permission.                                    #
#                                                                  #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           #
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      #
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         #
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         #
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            #
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     #
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         #
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; #
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER #
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      #
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      #
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       #
####################################################################


PROJECT(runUnitTests_singlecpu)

FIND_PACKAGE(Threads REQUIRED)

LIST(APPEND READDY_SINGLECPU_TEST_SOURCES TestMain.cpp)
LIST(APPEND READDY_SINGLECPU_TEST_SOURCES TestReactions.cpp)
LIST(APPEND READDY_SINGLECPU_TEST_SOURCES TestNeighborList.cpp)

ADD_EXECUTABLE(${PROJECT_NAME} ${READDY_SINGLECPU_TEST_SOURCES} ${GOOGLETEST_INCLUDE} ${READDY_INCLUDE_DIRS} ${TESTING_INCLUDE_DIR})
TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PUBLIC ${READDY_INCLUDE_DIRS} ${GOOGLETEST_INCLUDE} ${GOOGLEMOCK_INCLUDE} ${gtest_SOURCE_DIR} ${TESTING_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} readdy gtest_main Threads::Threads)
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES
        LINK_FLAGS "${EXTRA_LINK_FLAGS}"
        COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}")

ADD_TEST(NAME readdy-tests-singlecpu COMMAND ${PROJECT_NAME})

IF (READDY_INSTALL_UNIT_TEST_EXECUTABLE)
    INSTALL(TARGETS runUnitTests_singlecpu
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
ENDIF (READDY_INSTALL_UNIT_TEST_EXECUTABLE)