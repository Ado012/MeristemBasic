TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        src/common/myConfig.cc \
        src/common/myFiles.cc \
        src/common/myRandom.cc \
        src/common/mySignal.cc \
        src/common/myTimes.cc \
        src/compartment.cc \
        src/compartment/baseCompartmentNeighborhood.cc \
        src/compartment/compartmentNeighborhood.cc \
        src/cost/baseCostFunction.cc \
        src/cost/costFunction.cc \
        src/main/simulator.cc \
        src/nucleus.cc \
        src/organism.cc \
        src/reactions/baseReaction.cc \
        src/reactions/extendedMeristemReactions.cc \
        src/reactions/extendedMeristemReactionsHelperFunctions.cc \
        src/reactions/massAction.cc \
        src/solvers/baseSolver.cc \
        src/solvers/rungeKutta.cc \
        src/species.cc \
        src/topology.cc

HEADERS += \
    src/common/myConfig.h \
    src/common/myFiles.h \
    src/common/myRandom.h \
    src/common/mySignal.h \
    src/common/myTimes.h \
    src/common/typedefs.h \
    src/compartment.h \
    src/compartment/baseCompartmentNeighborhood.h \
    src/compartment/compartmentNeighborhood.h \
    src/cost/baseCostFunction.h \
    src/cost/costFunction.h \
    src/nucleus.h \
    src/organism.h \
    src/reactions/baseReaction.h \
    src/reactions/extendedMeristemReactions.h \
    src/reactions/extendedMeristemReactionsHelperFunctions.h \
    src/reactions/massAction.h \
    src/solvers/baseSolver.h \
    src/solvers/baseSolver.h \
    src/solvers/rungeKutta.h \
    src/solvers/rungeKutta.h \
    src/species.h \
    src/topology.h
