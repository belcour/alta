TEMPLATE        = lib
CONFIG         *=  plugin \
						 eigen

DESTDIR         = ../../build

INCLUDEPATH    += ../..
HEADERS         = data.h
SOURCES         = data.cpp


LIBS           += -L../../build   \
						-lcore          \
						-lfreeimageplus

