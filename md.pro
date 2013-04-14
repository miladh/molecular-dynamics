TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo -lconfig++


TARGET = md

release {
    DEFINES += ARMA_NO_DEBUG
    QMAKE_LFLAGS -= -O1
    QMAKE_LFLAGS += -O3
    QMAKE_LFLAGS_RELEASE -= -O1
    QMAKE_LFLAGS_RELEASE += -O3
    QMAKE_CXXFLAGS -= -O2
    QMAKE_CXXFLAGS += -O3
    QMAKE_CXXFLAGS_RELEASE -= -O2
    QMAKE_CXXFLAGS_RELEASE += -O3
}

SOURCES += src/main.cpp \
    src/fileManager/filemanager.cpp \
    src/System/system.cpp \
    src/atom/atom.cpp \
    src/generator/generator.cpp \
    src/force/twobodyforce.cpp \
    src/force/lj.cpp \
    src/force/noforce.cpp \
    src/modifier/modifier.cpp \
    src/modifier/andersenthermostat.cpp \
    src/modifier/berendsenthermostat.cpp \
    src/mdApp/mdapp.cpp\
    src/pores/pores.cpp \
    src/pores/circularpores.cpp \
    src/pores/cylindricalpores.cpp \
    src/force/onebodyforce.cpp \
    src/force/constantforce.cpp \
    src/force/force.cpp \
    src/changeDensity/changedensity.cpp \
    src/analysis/analysis.cpp

HEADERS += \
    src/includes/defines.h \
    src/fileManager/filemanager.h \
    src/System/system.h \
    src/atom/atom.h \
    src/generator/generator.h \
    src/force/twobodyforce.h \
    src/force/lj.h \
    src/force/noforce.h \
    src/modifier/modifier.h \
    src/modifier/andersenthermostat.h \
    src/modifier/berendsenthermostat.h \
    src/mdApp/mdapp.h\
    src/pores/pores.h \
    src/pores/circularpores.h \
    src/pores/cylindricalpores.h \
    src/force/onebodyforce.h \
    src/force/constantforce.h \
    src/force/force.h \
    src/changeDensity/changedensity.h \
    src/analysis/analysis.h




OTHER_FILES += \
    src/config.cfg

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

#C++11 features
COMMON_CXXFLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS

