/****************************************************************************
** Meta object code from reading C++ file 'data_def.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../data_def.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#include <QtCore/QVector>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'data_def.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_DebyeFitThread_t {
    QByteArrayData data[11];
    char stringdata0[177];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_DebyeFitThread_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_DebyeFitThread_t qt_meta_stringdata_DebyeFitThread = {
    {
QT_MOC_LITERAL(0, 0, 14), // "DebyeFitThread"
QT_MOC_LITERAL(1, 15, 15), // "sig_writelogext"
QT_MOC_LITERAL(2, 31, 0), // ""
QT_MOC_LITERAL(3, 32, 11), // "std::string"
QT_MOC_LITERAL(4, 44, 22), // "sig_fittingthread_done"
QT_MOC_LITERAL(5, 67, 28), // "sig_plot_current_thread_Data"
QT_MOC_LITERAL(6, 96, 15), // "QVector<double>"
QT_MOC_LITERAL(7, 112, 20), // "sig_update_statusbar"
QT_MOC_LITERAL(8, 133, 21), // "sig_live_update_model"
QT_MOC_LITERAL(9, 155, 4), // "stop"
QT_MOC_LITERAL(10, 160, 16) // "receive_ui_ready"

    },
    "DebyeFitThread\0sig_writelogext\0\0"
    "std::string\0sig_fittingthread_done\0"
    "sig_plot_current_thread_Data\0"
    "QVector<double>\0sig_update_statusbar\0"
    "sig_live_update_model\0stop\0receive_ui_ready"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_DebyeFitThread[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       5,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   49,    2, 0x06 /* Public */,
       4,    1,   52,    2, 0x06 /* Public */,
       5,    1,   55,    2, 0x06 /* Public */,
       7,    1,   58,    2, 0x06 /* Public */,
       8,    3,   61,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       9,    0,   68,    2, 0x0a /* Public */,
      10,    0,   69,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, 0x80000000 | 3,    2,
    QMetaType::Void, QMetaType::Int,    2,
    QMetaType::Void, 0x80000000 | 6,    2,
    QMetaType::Void, QMetaType::Double,    2,
    QMetaType::Void, 0x80000000 | 6, 0x80000000 | 6, 0x80000000 | 6,    2,    2,    2,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void DebyeFitThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        DebyeFitThread *_t = static_cast<DebyeFitThread *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->sig_writelogext((*reinterpret_cast< std::string(*)>(_a[1]))); break;
        case 1: _t->sig_fittingthread_done((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: _t->sig_plot_current_thread_Data((*reinterpret_cast< QVector<double>(*)>(_a[1]))); break;
        case 3: _t->sig_update_statusbar((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 4: _t->sig_live_update_model((*reinterpret_cast< QVector<double>(*)>(_a[1])),(*reinterpret_cast< QVector<double>(*)>(_a[2])),(*reinterpret_cast< QVector<double>(*)>(_a[3]))); break;
        case 5: _t->stop(); break;
        case 6: _t->receive_ui_ready(); break;
        default: ;
        }
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        switch (_id) {
        default: *reinterpret_cast<int*>(_a[0]) = -1; break;
        case 0:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< std::string >(); break;
            }
            break;
        case 2:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QVector<double> >(); break;
            }
            break;
        case 4:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 2:
            case 1:
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QVector<double> >(); break;
            }
            break;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (DebyeFitThread::*_t)(std::string );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&DebyeFitThread::sig_writelogext)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (DebyeFitThread::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&DebyeFitThread::sig_fittingthread_done)) {
                *result = 1;
                return;
            }
        }
        {
            typedef void (DebyeFitThread::*_t)(QVector<double> );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&DebyeFitThread::sig_plot_current_thread_Data)) {
                *result = 2;
                return;
            }
        }
        {
            typedef void (DebyeFitThread::*_t)(double );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&DebyeFitThread::sig_update_statusbar)) {
                *result = 3;
                return;
            }
        }
        {
            typedef void (DebyeFitThread::*_t)(QVector<double> , QVector<double> , QVector<double> );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&DebyeFitThread::sig_live_update_model)) {
                *result = 4;
                return;
            }
        }
    }
}

const QMetaObject DebyeFitThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_DebyeFitThread.data,
      qt_meta_data_DebyeFitThread,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *DebyeFitThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *DebyeFitThread::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_DebyeFitThread.stringdata0))
        return static_cast<void*>(const_cast< DebyeFitThread*>(this));
    return QThread::qt_metacast(_clname);
}

int DebyeFitThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    }
    return _id;
}

// SIGNAL 0
void DebyeFitThread::sig_writelogext(std::string _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void DebyeFitThread::sig_fittingthread_done(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void DebyeFitThread::sig_plot_current_thread_Data(QVector<double> _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void DebyeFitThread::sig_update_statusbar(double _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void DebyeFitThread::sig_live_update_model(QVector<double> _t1, QVector<double> _t2, QVector<double> _t3)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)), const_cast<void*>(reinterpret_cast<const void*>(&_t3)) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}
QT_END_MOC_NAMESPACE
