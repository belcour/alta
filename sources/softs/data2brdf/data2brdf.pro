CONFIG         += console \
						eigen
CONFIG         -= app_bundle

DESTDIR         = ../../build
INCLUDEPATH    += ../../ ../../libs/rational_1d \

SOURCES        += main.cpp
LIBS           += -L../../build -lcore

unix{
   PRE_TARGETDEPS += ../../build/libcore.a
	LIBS += -ldl
}
