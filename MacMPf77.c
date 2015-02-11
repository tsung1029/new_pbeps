/* This is a Fortran callable interface library to the C version of a
   Multitasking library for Apple Macintosh
   based on the Multiprocessing Library in the MacOS, v. 2.1
   The file MPerrs is used throughout for error messages
   written by viktor k. decyk, ucla
   copyright 2000, regents of the university of california.
   all rights reserved.
   no warranty for proper operation of this software is given or implied.
   software or information may be copied, distributed, and used at own
   risk; it may not be distributed without this notice included verbatim
   with each file.  Designed for use with 64 bit addressing.
   update: August 6, 2009                                             */

#include "MacMP.h"

void MP_Cr8task(void *(*entryPoint) (void *), long *ttype, void (*proc)(),
                void *parameter, long *nargs, long *taskid);

void mp_init_(long *nproc) {
   MP_Init(nproc);
   return;
}

void mp_cr8task_(long (*entryPoint)(void *parameter), long *ttype,
                 void (*proc)(), void *parameter, long *nargs,
                 long *taskid) {
   MP_Cr8task((void *)entryPoint,ttype,proc,parameter,nargs,taskid);
   return;
}

void mp_taskwait_(long *taskid) {
   MP_Taskwait(taskid);
   return;
}

int mp_sndsig_(long *taskid) {
   return MP_Sndsig(taskid);
}

int mp_waitsig_(long *taskid) {
   return MP_Waitsig(taskid);
}

void mp_killtask_(long *taskid) {
   MP_Killtask(taskid);
   return;
}

void mp_end_() {
   MP_End();
   return;
}

void mp_taskstart_(long *taskid, void (*proc)(), long *nargs, long *arg1,
                  long *arg2, long *arg3, long *arg4, long *arg5, long *arg6,
                  long *arg7, long *arg8, long *arg9, long *arg10, long *arg11,
                  long *arg12, long *arg13, long *arg14, long *arg15,
                  long *arg16, long *arg17, long *arg18, long *arg19,
                  long *arg20, long *arg21, long *arg22, long *arg23,
                  long *arg24, long *arg25) {
   switch (*nargs) {
   case 0:
      MP_Taskstart(taskid,proc,nargs);
      return;
   case 1:
      MP_Taskstart(taskid,proc,nargs,arg1);   
      return;
   case 2:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2);
      return;
   case 3:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3);
      return;
   case 4:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4);
      return;
   case 5:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5);
      return;
   case 6:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6);
      return;
   case 7:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7);
      return;
   case 8:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8);
      return;
   case 9:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9);
      return;
   case 10:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10);
      return;
   case 11:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11);
      return;
   case 12:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12);
      return;
   case 13:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13);
      return;
   case 14:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14);
      return;
   case 15:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15);
      return;
   case 16:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16);
      return;
   case 17:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17);
      return;
   case 18:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18);
      return;
   case 19:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19);
      return;
   case 20:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20);
      return;
   case 21:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21);
      return;
   case 22:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22);
      return;
   case 23:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22,arg23);
      return;
   case 24:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24);
      return;
   case 25:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24,
                   arg25);
      return;
   }
   return;
}

void mp_taskinit_(long *taskid, void (*proc)(), long *nargs, long *arg1,
                  long *arg2, long *arg3, long *arg4, long *arg5, long *arg6,
                  long *arg7, long *arg8, long *arg9, long *arg10, long *arg11,
                  long *arg12, long *arg13, long *arg14, long *arg15,
                  long *arg16, long *arg17, long *arg18, long *arg19,
                  long *arg20, long *arg21, long *arg22, long *arg23,
                  long *arg24, long *arg25) {

   switch (*nargs) {
   case 0:
      MP_Taskinit(taskid,proc,nargs);
      return;
   case 1:
      MP_Taskinit(taskid,proc,nargs,arg1);   
      return;
   case 2:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2);
      return;
   case 3:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3);
      return;
   case 4:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4);
      return;
   case 5:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5);
      return;
   case 6:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6);
      return;
   case 7:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7);
      return;
   case 8:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8);
      return;
   case 9:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9);
      return;
   case 10:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10);
      return;
   case 11:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11);
      return;
   case 12:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12);
      return;
   case 13:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13);
      return;
   case 14:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14);
      return;
   case 15:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15);
      return;
   case 16:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16);
      return;
   case 17:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17);
      return;
   case 18:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18);
      return;
   case 19:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19);
      return;
   case 20:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20);
      return;
   case 21:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21);
      return;
   case 22:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22);
      return;
   case 23:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22,arg23);
      return;
   case 24:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24);
      return;
   case 25:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24,
                   arg25);
      return;
   }
   return;
}

void mp_setstack_(long *stackval) {
   long lstackval;
   lstackval = *stackval;
   MP_Setstack(lstackval);
   return;
}

void mp_taskbuild_(long *taskid, void (*proc)(), long *nargs, long *arg1,
                  long *arg2, long *arg3, long *arg4, long *arg5, long *arg6,
                  long *arg7, long *arg8, long *arg9, long *arg10, long *arg11,
                  long *arg12, long *arg13, long *arg14, long *arg15,
                  long *arg16, long *arg17, long *arg18, long *arg19,
                  long *arg20, long *arg21, long *arg22, long *arg23,
                  long *arg24, long *arg25) {
   switch (*nargs) {
   case 0:
      MP_Taskbuild(taskid,proc,nargs);
      return;
   case 1:
      MP_Taskbuild(taskid,proc,nargs,arg1);   
      return;
   case 2:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2);
      return;
   case 3:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3);
      return;
   case 4:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4);
      return;
   case 5:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5);
      return;
   case 6:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6);
      return;
   case 7:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7);
      return;
   case 8:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8);
      return;
   case 9:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9);
      return;
   case 10:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10);
      return;
   case 11:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11);
      return;
   case 12:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12);
      return;
   case 13:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13);
      return;
   case 14:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14);
      return;
   case 15:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15);
      return;
   case 16:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16);
      return;
   case 17:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17);
      return;
   case 18:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18);
      return;
   case 19:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19);
      return;
   case 20:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20);
      return;
   case 21:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21);
      return;
   case 22:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22);
      return;
   case 23:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22,arg23);
      return;
   case 24:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24);
      return;
   case 25:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24,
                   arg25);
      return;
   }
   return;
}

void mp_runtask_(long *taskid) {
   MP_Runtask(taskid);
   return;
}

void mp_initialized_(long *flag) {
   MP_Initialized(flag);
   return;
}
