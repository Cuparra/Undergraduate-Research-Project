TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    ParticleSwarmOptimization.c \
    ParetoArchive.c \
    main.c \
    Function.c

HEADERS += \
    ParetoArchive.h \
    Swarm.h \
    ParticleSwarmOptimization.h \
    Function.h

