######################################################################
# Automatically generated by qmake (3.1) Thu Oct 7 16:52:33 2021
######################################################################

QT += widgets printsupport

TEMPLATE = app
TARGET = ksemawc

# You can make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# Please consult the documentation of the deprecated API in order to know
# how to port your code away from it.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

win32 {
    QT += core gui
    include( .\QWT_win.pri )
    CONFIG+= c++11 gui console
    LIBS += -L"C:/Program Files/Tools/mingw810_64/cminpack-1.3.8/build" -lcminpack
    HEADERS += cminpack.h
}

unix {
    CONFIG += link_pkgconfig
    PKGCONFIG += cminpack cblas blas qwt
}
INCLUDEPATH += .

# Input
HEADERS += ksemawc.h
FORMS += ksemawc.ui
SOURCES += ksemawc.cpp main.cpp

target.path=.\
INSTALLS += target

