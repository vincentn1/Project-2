TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
INCLUDEPATH += C:\Users\Johannes\C++\Libraries\armadillo-5.500.2\include
QMAKE_CXXFLAGS += -O3

LIBS += -L C:\Users\Johannes\C++\Libraries\lapack_blas\32bit

LIBS += -llibblas -lliblapack

SOURCES += main.cpp

