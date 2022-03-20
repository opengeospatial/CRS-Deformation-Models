/***************************************************************************
 $Id$

 Copyright 2011 Crown copyright (c)
 Land Information New Zealand and the New Zealand Government.
 All rights reserved

 This program is released under the terms of the new BSD license. See 
 the LICENSE file for more information.
****************************************************************************/

#include <memory.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "get_image_path.h"

#if !defined(_WIN32) && !defined(_MSC_VER)
#include <sys/types.h>
#include <unistd.h>

char* get_image_path()
{
    char* path;
    char _link[20];
    char buf[10];
    pid_t pid = getpid();
    sprintf( buf,"%d", pid );
    strcpy( _link, "/proc/" );
    strcat( _link, buf );
#if defined(__linux) || defined(linux)
    strcat( _link, "/exe" );
#endif
#if defined(sun) || defined(__sun)
    strcat( _link, "/path/a.out" );
#endif
#if defined(__bsdi__)
    strcat( _link, "/file" );
#endif
    char proc[512];
    ssize_t len = readlink( _link, proc, 512);
    if ( len != -1 )
    {
      proc[len] = '\0';
      path = new char[strlen( proc ) + 1];
      path = strcpy( path, proc );
    }
    return path;
}
#else

char *get_image_path()
{
    const char *path = _pgmptr;
    int len = strlen(path);
    if( len > 4 && _stricmp(path+len-4,".exe") == 0 )
    {
        len -= 4;
    }
    char *image_path = new char[len+1];
    strncpy(image_path, path, len );
    image_path[len] = 0;
    return image_path;

}

#endif
