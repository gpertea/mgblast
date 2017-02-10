/*
 * $Id: MW-BSD-limits-fix.h,v 1.1 2005/08/31 19:03:13 rsmith Exp $
 *
 * CodeWarrior prefix file for BSD builds to fix bug in Apple's machine/limits.h
 */
#ifdef __MWERKS__
# ifndef __CHAR_BIT__ 
#  define __CHAR_BIT__ 8
# endif
# ifndef __SCHAR_MAX__ 
#  define __SCHAR_MAX__ 255
# endif
#endif