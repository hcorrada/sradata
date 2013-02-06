AC_INIT(DESCRIPTION)
AC_PROG_CC

AC_ARG_WITH([sra-includes], 
  [AS_HELP_STRING([--with-sra-includes],[Use SRA includes in DIR])],
  SRAINCLUDES=$withval,
  SRAINCLUDES="")

AC_ARG_WITH([sra-libs],
  [AS_HELP_STRING([--with-sra-libs],[USE SRA libraries in DIR])],
  SRALIBS=$withval,
  SRALIBS="")

if test [ -z "$SRAINCLUDES" -or -z "$SRALIBS" ] ; then
   AC_MSG_ERROR([Locations of SRA headers and libraries MUST be given])
fi

if test [ ! -d "$SRAINCLUDES" -or ! -d "$SRALIBS" ] ; then
   AC_MSG_ERROR([Locations of SRA headers and libraries not valid])
fi

PKG_CPPFLAGS="-I${SRAINCLUDES}"
CPPFLAGS=${PKG_CPPFLAGS}
AC_CHECK_HEADER(sra/sradb.h)
if test "${ac_cv_header_sra_sradb_h}" == no ; then
   AC_MSG_ERROR([SRA header sra/sradb.h not found in ${SRAINCLUDES}])
fi

PKG_LIBS="-L${SRALIBS} -lsradb -lvdb -lklib -lkascii -lm -lz -lbz2 -ldl -lpthread"

AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_LIBS)
AC_OUTPUT(src/Makevars)


