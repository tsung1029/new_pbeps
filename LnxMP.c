/* Multitasking library based on Pthreads in Linux
   The file MPerrs is used throughout for error messages
   written by viktor k. decyk, ucla
   copyright 2003, regents of the university of california.
   all rights reserved.
   no warranty for proper operation of this software is given or implied.
   software or information may be copied, distributed, and used at own
   risk; it may not be distributed without this notice included verbatim
   with each file. Designed for use with 64 bit addressing.
   update: August 6, 2009                                               */

#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include "MacMP.h"

#include <sys/types.h>
#include <pthread.h>

/* MAXARGS = maximum number of arguments in task procedure */
#define MAXARGS                 25
/* MAXTASKS = maximum number of tasks supported */
#define MAXTASKS                16

/* internal common block for tasks
   notifyq = (0,1) = (no,yes) multi-tasking initialized
   ntask = task record
   itask = pthread_t record
   params = arguments in task procedures
   smutex = mutex variable for condition variables
   mqueues = condition variables for tasks (1=request,2=result)
   mvalues = semaphore variables for tasks (1=request,2=result)
   dfsz = default stack size
   stackSize = size of stack for task in bytes (0 = default)  */

   static long notifyq = 0;
   static int ntask[MAXTASKS];
   static pthread_t itask[MAXTASKS];
   static void *params[MAXTASKS][MAXARGS+6];
   static pthread_mutex_t smutex;
   static pthread_cond_t mqueues[MAXTASKS][2];
   static int mvalues[MAXTASKS][2];
   static size_t dfsz = 0;
   static pthread_attr_t stackSize;

static FILE *unit2;

/* prototypes for internal procedures */

void *tproc(void *parameter);

void *tsproc(void *parameter);

void  mproc(void (*proc)(), long nargs, void *parameter);

int MPProcessors();

void MP_Init(long *nproc) {
/* initialize multitasking environment
   nproc = number of processors, 0 if Multitasking is not available
   output: nproc
local data                             */
   int i, oss;
/* internal common block for tasks
   notifyq = (0,1) = (no,yes) multi-tasking initialized
   ntask = task record
   itask = pthread_t record
   smutex = mutex variable for condition variables
   mvalues = semaphore variables for tasks (1=request,2=result)
   dfsz = default stack size
   stackSize = size of stack for task in bytes (0 = default)    */
   *nproc = 0;
/* Open error file */
   unit2 = fopen("MPerrs","w");
/* Return the number of processors on the host computer */
   *nproc = MPProcessors();
/* *nproc = LPProcessors(); */
/* Clear task records */
   for (i = 0; i < MAXTASKS; i++) {
      ntask[i] = 0;
      itask[i] = 0;
      mvalues[i][0] = -1;
      mvalues[i][1] = -1;
   }
/* Initialize mutex */
   oss = pthread_mutex_init(&smutex,NULL);
   if (oss)
      fprintf(unit2,"pthread_mutex_init error = %d,%s\n",oss,
              strerror(oss));
/* set current stacksize (in bytes) for tasks */
   oss = pthread_attr_init(&stackSize);
   if (oss)
      fprintf(unit2,"pthread_attr_init error = %d,%s\n",oss,
              strerror(oss));
   oss = pthread_attr_getstacksize(&stackSize,&dfsz);
   if (oss)
      fprintf(unit2,"pthread_attr_getstacksize error = %d,%s\n",
              oss,strerror(oss));
   oss = pthread_attr_setstacksize(&stackSize,(size_t)16384);
   if (oss)
      fprintf(unit2,"pthread_attr_setstacksize error = %d,%s\n",
              oss,strerror(oss));
/* MP Initialized */
   notifyq = 1;
   return;
}

void MP_Cr8task(void *(*entryPoint) (void *), long *ttype, void (*proc)(),
               void *parameter, long *nargs, long *taskid) {
/* start task
   overhead is about 130 microseconds on Macintosh G4/450
   entrypoint = pointer to task procedure
   ttype = task type (0=non-reusable,1=reusable(signaling),2=non-running)
   proc = pointer to actual procedure
   parameter = pointer to task procedure argument
   nargs = number of arguments in task procedure
   taskid = index to task record (0 if error)
   input: all except taskid
   output: taskid
local data                                              */
   void **param;
   int i, j, oss;
/* internal common block for tasks
   notifyq = (0,1) = (no,yes) multi-tasking initialized
   ntask = task record
   itask = pthread_t record
   params = arguments in task procedures
   mqueues = condition variables for tasks (1=request,2=result)
   mvalues = semaphore variables for tasks (1=request,2=result)
   stackSize = size of stack for task in bytes (0 = default)  */
   param = (void **)parameter;
   *taskid = 0;
/* Check for errors */
   if ((!notifyq) && (*ttype < 2)) {
      fprintf(unit2,"MP_Cr8task: MP not initialized\n");
      return;
   }
   else if (*nargs > MAXARGS) {
      fprintf(unit2,"Too many arguments in task, nargs = %ld, MAXARGS = %d\n",
              *nargs,MAXARGS);
      return;
   }
/* Find space for record */
   i = -1;
L10: i += 1;
   if (i >= MAXTASKS) {
      fprintf(unit2,"Exceeded maximum number of tasks = %d\n",MAXTASKS);
      return;
   }
   else if (ntask[i])
      goto L10;
/* Create semaphores */
   if (*ttype==1) {
      for (j = 0; j < 2; j++) {
         oss = pthread_cond_init(&mqueues[i][j],NULL);
         if (oss) {
            fprintf(unit2,"pthread_cond_init error = %d,%s\n",oss,
                 strerror(oss));
            return;
         }
         mvalues[i][j] = 0;
      }
   }
   ntask[i] = i + 1;
/* Copy number of arguments */
   params[i][0] = (void *)(*nargs);
/* Copy arguments */
   for (j = 0; j < *nargs; j++) {
      params[i][j+1] = param[j];
   }
/* Copy procedure name and message queue pointers */
   params[i][*nargs+1] = (void *)proc;
   params[i][*nargs+2] = (void *)&mqueues[i][0];
   params[i][*nargs+3] = (void *)&mqueues[i][1];
   params[i][*nargs+4] = (void *)&mvalues[i][0];
   params[i][*nargs+5] = (void *)&mvalues[i][1];
/* Create a preemptive task */
   if (*ttype < 2) {
      oss = pthread_create(&itask[i],&stackSize,entryPoint,
                           &params[i][0]);
      if (oss) {
         fprintf(unit2,"pthread_create error = %d, %s\n",oss,
                 strerror(oss));
         ntask[i] = 0;
         itask[i] = 0;
         return;
      }
   }
   *taskid = i + 1;
   return;
}

void MP_Taskwait(long *taskid) {
/* wait for task to complete
   taskid = index to task record (0 if task completed successfully)
   input and output: taskid
local data                                             */
   int i, j, oss;
   void *result = 0;
/* internal common block for tasks
   notifyq = (0,1) = (no,yes) multi-tasking initialized
   ntask = task record
   itask = pthread_t record
   mqueues = condition variables for tasks (1=request,2=result)
   mvalues = semaphore variables for tasks (1=request,2=result) */
   i = *taskid - 1;
/* Check for errors */
   if (!notifyq) {
      fprintf(unit2,"MP_Taskwait: MP not initialized\n");
      return;
   }
   else if ((i < 0) || (i >=MAXTASKS)) {
      fprintf(unit2,"MP_Taskwait: Invalid taskid = %ld\n",*taskid);
      return;
   }
/* Look for message from taskid */
   if (ntask[i] > 0) {
/* Wait for specified thread's termination */
      oss = pthread_join(itask[i],&result);
      if (oss) {
         fprintf(unit2,"pthread_join error, taskid, oss = %d,%d,%s\n",
                 i+1,oss,strerror(oss));
         return;
      }
/* Mark task as completed */
      ntask[i] = 0;
   }
/* Remove semaphores */
   for (j = 0; j < 2; j++) {
      if (mvalues[i][j] >= 0) {
         oss = pthread_cond_destroy(&mqueues[i][j]);
         if (oss)
            fprintf(unit2,"pthread_cond_destroy error, taskid,oss=%d,%d,%s\n",
                    i+1,oss,strerror(oss));
         mvalues[i][j] = -1;
      }
   }
   itask[i] = 0;
   *taskid = 0;
   return;
}

int MP_Sndsig(long *taskid) {
/* send message to task
   overhead for sending a message to a task and receiving a reply 
   is about 60 microseconds on Macintosh G4/450
   taskid = index to task record
   input: all
local data                                          */
   int i, oss, err;
/* internal common block for tasks
   notifyq = (0,1) = (no,yes) multi-tasking initialized
   ntask = task record
   mqueues = condition variables for tasks (1=request,2=result)
   mvalues = semaphore variables for tasks (1=request,2=result) */
   i = *taskid - 1;
/* Check for errors */
   if (!notifyq) {
      fprintf(unit2,"MP_Sndsig: MP not initialized\n");
      return 1;
   }
   else if ((i < 0) || (i >=MAXTASKS)) {
      fprintf(unit2,"MP_Sndsig: Invalid taskid = %ld\n",*taskid);
      return 2;
   }
   else if (!(ntask[i])) {
      fprintf(unit2,"MP_Sndsig: Task already ended, taskid = %ld\n",
              *taskid);
      return 3;
   }
   else if (mvalues[i][0] < 0) {
      fprintf(unit2,"MP_Sndsig: Invalid type for task = %ld\n",*taskid);
      return 4;
   }
/* Signal a semaphore */
   err = pthread_mutex_lock(&smutex);
   mvalues[i][0] = 1;
   oss = pthread_cond_signal(&mqueues[i][0]);
   err = pthread_mutex_unlock(&smutex);
   if (oss) {
      fprintf(unit2,"pthread_cond_signal error for taskid,oss=%d,%d,%s\n",
              i+1,oss,strerror(oss));
      return 5;
   }
   return 0;
}

int MP_Waitsig(long *taskid) {
/* receive signal from task
   overhead for sending a message to a task and receiving a reply 
   is about 60 microseconds on Macintosh G4/450
   taskid = index to task record
   input: taskid
   input: all
local data                                          */
   int i, oss, err;
/* internal common block for tasks
   notifyq = (0,1) = (no,yes) multi-tasking initialized
   ntask = task record
   mqueues = condition variables for tasks (1=request,2=result) */
   i = *taskid - 1;
/* Check for errors */
   if (!notifyq) {
      fprintf(unit2,"MP_Waitsig: MP not initialized\n");
      return 1;
   }
   else if ((i < 0) || (i >=MAXTASKS)) {
      fprintf(unit2,"MP_Waitsig: Invalid taskid = %ld\n",*taskid);
      return 2;
   }
   else if (!(ntask[i])) {
      fprintf(unit2,"MP_Waitsig: Task already ended, taskid = %ld\n",
              *taskid);
      return 3;
   }
   else if (mvalues[i][1] < 0) {
      fprintf(unit2,"MP_Waitsig: Invalid type for task = %ld\n",*taskid);
      return 4;
   }
/* Wait on semaphore */
   err = pthread_mutex_lock(&smutex);
   if (mvalues[i][1]==0) {
      oss = pthread_cond_wait(&mqueues[i][1],&smutex);
      if (oss) {
         fprintf(unit2,"pthread_cond_wait error for taskid,oss=%d,%d,%s\n",
                 i+1,oss,strerror(oss));
         return 5;
      }
   }
   mvalues[i][1] = 0;
   err = pthread_mutex_unlock(&smutex);
   return 0;
}

void MP_Killtask(long *taskid) {
/* terminate a task
   taskid = index to task record (0 if task killed successfully)
local data                             */
   int i, oss;
/* internal common block for tasks
   notifyq = (0,1) = (no,yes) multi-tasking initialized
   itask = pthread_t record                */
   i = *taskid - 1;
/* Check for errors */
   if (!notifyq) {
      fprintf(unit2,"MP_Killtask: MP not initialized\n");
      return;
   }
   else if ((i < 0) || (i >=MAXTASKS)) {
      fprintf(unit2,"MP_Killtask: Invalid taskid = %ld\n",*taskid);
      return;
   }
   else if (!(ntask[i])) {
      fprintf(unit2,"Task already killed, taskid = %ld\n",*taskid);
      return;
   }
/* Terminate an existing task */
   oss = pthread_cancel(itask[i]);
   if (oss)
      fprintf(unit2,"pthread_cancel error, taskid, oss = %d,%d,%s\n",
              i+1,oss,strerror(oss));
   else {
 /* Wake up reusable task */
      if (mvalues[i][0] >= 0)
         oss = MP_Sndsig(taskid);
      MP_Taskwait(taskid);
   }
   return;
}

void MP_End() {
/* terminate multitasking environment
local data                             */
   int i;
   long taskid;
   long oss;
/* internal common block for tasks
   notifyq = notifyqueue record for tasks */
/* Check for errors */
   if (!notifyq) {
      fprintf(unit2,"MP_End: MP not initialized\n");
      return;
   }
/* Check if any tasks are outstanding */
   for (i = 0; i < MAXTASKS; i++) {
      if (ntask[i] > 0) {
         taskid = i + 1;
         fprintf(unit2,"MP_End: Task still outstanding, taskid=%d\n",i+1);
         MP_Killtask(&taskid);
      }
   }
/* Destroy mutex */
   oss = pthread_mutex_destroy(&smutex);
   if (oss)
      fprintf(unit2,"pthread_mutex_destroy error = %ld,%s\n",oss,
              strerror(oss));
/* Destroy attribute */
   oss = pthread_attr_destroy(&stackSize);
   if (oss)
      fprintf(unit2,"pthread_attr_destroy error = %ld,%s\n",oss,
              strerror(oss));
/* Delete file if empty */
   if (!fseek(unit2,0,SEEK_END)) {
      oss = ftell(unit2);
      fclose(unit2);
      if (!oss)
         remove("MPerrs");
   }
   notifyq = 0;
   return;
}

void MP_Taskstart(long *taskid, void (*proc)(), long *nargs, ...) {
/* create a task by packing arguments into a single structure
   taskid = index to notify queue for task (0 if error)
   proc = pointer to actual procedure
   nargs = number of arguments in actual procedure
   restrictions: actual arguments should be pointers and not temporary
   variables, such as return values of functions.  Also, local variables
   in procedure which called MP_Taskstart should not be used if
   MP_Taskwait is not called in the same procedure.
local data                                                           */
   void *param[MAXARGS];
   int i;
   long ttype = 0;
   va_list argptr;
   va_start(argptr,nargs);
/* Create argument list for task */
   for (i = 0; i < *nargs; i++) {
      param[i] = va_arg(argptr,void*);
   }
   va_end(argptr);
/* Start task */
   MP_Cr8task(&tproc,&ttype,proc,&param[0],nargs,taskid);
   if (!(*taskid))
      fprintf(unit2,"MP_Taskstart failed\n");
   return;
}

void MP_Taskinit(long *taskid, void (*proc)(), long *nargs, ...) {
/* create a reusable task by packing arguments into a single structure
   taskid = index to notify queue for task (0 if error)
   proc = pointer to actual procedure
   nargs = number of arguments in actual procedure
   restrictions: actual arguments should be pointers and not temporary
   variables, such as return values of functions.  Also, local variables
   in procedure which called MP_Taskinit should not be used if
   MP_Sndsig and MP_Waitsig are not called in the same procedure.
local data                                                           */
   void *param[MAXARGS];
   int i;
   long ttype = 1;
   va_list argptr;
   va_start(argptr,nargs);
/* Create argument list for task */
   for (i = 0; i < *nargs; i++) {
      param[i] = va_arg(argptr,void*);
   }
   va_end(argptr);
/* Start task */
   MP_Cr8task(&tsproc,&ttype,proc,&param[0],nargs,taskid);
   if (!(*taskid))
      fprintf(unit2,"MP_Taskinit failed\n");
   return;
}

void MP_Setstack(long stackval) {
/* set current stacksize (in bytes) for tasks 
   input: stackval
local data                                    */
   int oss = 0;
/* internal common block for tasks
   dfsz = default stack size
   stackSize = size of stack for task in bytes (0 = default) */
   if (stackval==0)
      oss = pthread_attr_setstacksize(&stackSize,dfsz);
   else if (stackval > 0)
      oss = pthread_attr_setstacksize(&stackSize,(size_t)stackval);
   if (oss)
      fprintf(unit2,"pthread_attr_setstacksize, error = %ld,%d,%s\n",
              stackval,oss,strerror(oss));
   return;
}

void MP_Taskbuild(long *taskid, void (*proc)(), long *nargs, ...) {
/* create a non-running task by packing arguments into a single structure
   generally used for debugging in conjunction with MP_Runtask
   taskid = index to notify queue for task (0 if error)
   proc = pointer to actual procedure
   nargs = number of arguments in actual procedure
   restrictions: actual arguments should be pointers and not temporary
   variables, such as return values of functions.  Also, local variables
   in procedure which called MP_Taskbuild should not be used if
   MP_Runtask is not called in the same procedure.
local data                                                           */
   void *param[MAXARGS];
   int i;
   long ttype = 2;
   va_list argptr;
   va_start(argptr,nargs);
/* Create argument list for task */
   for (i = 0; i < *nargs; i++) {
      param[i] = va_arg(argptr,void*);
   }
   va_end(argptr);
/* Start task */
   MP_Cr8task(&tproc,&ttype,proc,&param[0],nargs,taskid);
   if (!(*taskid))
      fprintf(unit2,"MP_Taskbuild failed\n");
   return;
}

void MP_Runtask(long *taskid) {
/* executes a non-running task manually
   generally used for debugging in conjunction with MP_Taskbuild
   taskid = index to task record (0 if task completed successfully)
   input: taskid
local data  */
   int i, ierr;
/* internal common block for tasks
   params = arguments in task procedures */
   i = *taskid - 1;
/* Check for errors */
   if ((i < 0) || (i >=MAXTASKS)) {
      fprintf(unit2,"MP_Runtask: Invalid taskid = %ld\n",*taskid);
      return;
   }
   ierr = *(long *)tproc(&params[i][0]);
   ntask[i] = 0;
   *taskid = 0;
   return;
}

void MP_Initialized(long *flag) {
/* indicate whether MP_Init has been called
   flag = true if MP_Init has been called, false otherwise
   output: flag                                             */
/* internal common block for tasks
   notifyq = (0,1) = (no,yes) multi-tasking initialized     */
   if (notifyq)
      *flag = 1;
   else
      *flag = 0;
   return;
}

/*--------------------------------------------------------------------*/
/* Internal functions for Multitasking library                        */
/*--------------------------------------------------------------------*/

void mproc(void (*proc)(), long nargs, void *parameter);

void *tproc(void *parameter) {
/* generic task procedure
local data                */
   void **param;
   int otype, ostate, oss;
   long nargs;
/* Make thread cancelable */
   oss = pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS,&otype);
   oss = pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,&ostate);
   pthread_testcancel();
/* Extract arguments and execute */
   param = (void **)parameter;
   nargs = (long )param[0];
/* Run task */
   mproc(param[nargs+1],nargs,&param[1]);
   return (void *)0;
}

void *tsproc(void *parameter) {
/* generic task procedure which waits on semaphores
local data                */
   void **param;
   int otype, ostate, oss, err;
   long nargs;
/* Check if thread is canceled */
   oss = pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS,&otype);
   oss = pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,&ostate);
   pthread_testcancel();
/* Disable cancelation */
   oss = pthread_setcancelstate(PTHREAD_CANCEL_DISABLE,&ostate);
/* Extract arguments and execute */
   param = (void **)parameter;
   nargs = (long )param[0];
/* Wait on semaphore */
L10: err = pthread_mutex_lock(&smutex);
   if (*(int *)param[nargs+4]==0)
      oss = pthread_cond_wait(param[nargs+2],&smutex);
   else
      oss = 0;
   *(int *)param[nargs+4] = 0;
   err = pthread_mutex_unlock(&smutex);
/* Signal arrived */
   if (!oss) {
/* Check if thread is canceled */
      oss = pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,&ostate);
      pthread_testcancel();
/* Run task */
      mproc(param[nargs+1],nargs,&param[1]);
/* Disable cancelation */
      oss = pthread_setcancelstate(PTHREAD_CANCEL_DISABLE,&ostate);
/* Signal a semaphore */
      err = pthread_mutex_lock(&smutex);
      *(int *)param[nargs+5] = 1;
      oss = pthread_cond_signal(param[nargs+3]);
      err = pthread_mutex_unlock(&smutex);
/* Wait for next signal */
      if (!oss)
         goto L10;
   }
   return (void *)(long )oss;
}

void mproc(void (*proc)(), long nargs, void *parameter) {
/* task subroutine with multiple possible arguments
   proc = pointer to actual procedure
   nargs = number of arguments in task procedure
   parameter = pointer to task procedure argument */
   void **param;
   param = (void **)parameter;
   switch (nargs) {
      case 0:
         proc();
         break;
      case 1:
         proc((void *)param[0]);
         break;
      case 2:
         proc((void *)param[0],(void *)param[1]);
         break;
      case 3:
         proc((void *)param[0],(void *)param[1],(void *)param[2]);
         break;
      case 4:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3]);
         break;
      case 5:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4]);
         break;
      case 6:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5]);
         break;
      case 7:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6]);
         break;
      case 8:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7]);
         break;
      case 9:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8]);
         break;
      case 10:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9]);
         break;
      case 11:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10]);
         break;
      case 12:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11]);
         break;
      case 13:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12]);
         break;
      case 14:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13]);
         break;
      case 15:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14]);
         break;
      case 16:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15]);
         break;
      case 17:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16]);
         break;
      case 18:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16],(void *)param[17]);
         break;
      case 19:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16],(void *)param[17],
              (void *)param[18]);
         break;
      case 20:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16],(void *)param[17],
              (void *)param[18],(void *)param[19]);
         break;
      case 21:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16],(void *)param[17],
              (void *)param[18],(void *)param[19],(void *)param[20]);
         break;
      case 22:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16],(void *)param[17],
              (void *)param[18],(void *)param[19],(void *)param[20],
              (void *)param[21]);
         break;
      case 23:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16],(void *)param[17],
              (void *)param[18],(void *)param[19],(void *)param[20],
              (void *)param[21],(void *)param[22]);
         break;
      case 24:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16],(void *)param[17],
              (void *)param[18],(void *)param[19],(void *)param[20],
              (void *)param[21],(void *)param[22],(void *)param[23]);
      case 25:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16],(void *)param[17],
              (void *)param[18],(void *)param[19],(void *)param[20],
              (void *)param[21],(void *)param[22],(void *)param[23],
			  (void *)param[24]);
   }
   return;
}

void prparms(long taskid) {
/* debugging subroutine for printing task arguments
   used for debugging in conjunction with MP_Taskbuild
   taskid = index to task record (0 if task completed successfully)
   input: taskid
local data    */
   int i, n, nargs, iarg;
   void *larg;
   float arg;
   double darg;
/* internal common block for tasks
   params = arguments in task procedures */
   i = taskid - 1;
/* Check for errors */
   if ((i < 0) || (i >=MAXTASKS)) {
      return;
   }
/* extract number of arguments */
   nargs = (long )params[i][0];
   fprintf(unit2,"taskid, nargs = %ld,%d\n",taskid,nargs);
/* write location and values of arguments (as integer, real and double) */
   for (n = 0; n < nargs; n++) {
      larg = params[i][n+1];
      iarg = *(int *) larg;
      arg = *(float *) larg;
      darg = *(double *) larg;
      fprintf(unit2,"n, loc, int, float, double = %d,%p,%d,%f,%f\n",
	         n,larg,iarg,arg,darg);
   }
   return;
}
