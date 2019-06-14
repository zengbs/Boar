#include <stdio.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------*/
/*! \fn void* calloc_1d_array(size_t nc, size_t size)
 *  *  \brief Construct 1D array = array[nc]  */
void *
calloc_1d_array (size_t nc, size_t size)
{
  void *array;

  if ((array = (void *) calloc (nc, size)) == NULL)
    {
      printf ("[calloc_1d] failed to allocate memory (%d of size %d)\n", (int) nc, (int) size);
      return NULL;
    }
  return array;
}

/*----------------------------------------------------------------------------*/
/*! \fn void** calloc_2d_array(size_t nr, size_t nc, size_t size)
 *  *  \brief Construct 2D array = array[nr][nc]  */
void **
calloc_2d_array (size_t nr, size_t nc, size_t size)
{
  void **array;
  size_t i;

  if ((array = (void **) calloc (nr, sizeof (void *))) == NULL)
    {
      printf ("[calloc_2d] failed to allocate mem for %d pointers\n", (int) nr);
      return NULL;
    }

  if ((array[0] = (void *) calloc (nr * nc, size)) == NULL)
    {
      printf ("[calloc_2d] failed to allocate memory (%d X %d of size %d)\n", (int) nr, (int) nc, (int) size);
      free ((void *) array);
      return NULL;
    }

  for (i = 1; i < nr; i++)
    {
      array[i] = (void *) ((unsigned char *) array[0] + i * nc * size);
    }

  return array;
}


/*----------------------------------------------------------------------------*/
/*! \fn void*** calloc_3d_array(size_t nt, size_t nr, size_t nc, size_t size)
 *  *  \brief Construct 3D array = array[nt][nr][nc]  */
void ***
calloc_3d_array (size_t nt, size_t nr, size_t nc, size_t size)
{
  void ***array;
  size_t i, j;

  if ((array = (void ***) calloc (nt, sizeof (void **))) == NULL)
    {
      printf ("[calloc_3d] failed to allocate memory for %d 1st-pointers\n", (int) nt);
      return NULL;
    }

  if ((array[0] = (void **) calloc (nt * nr, sizeof (void *))) == NULL)
    {
      printf ("[calloc_3d] failed to allocate memory for %d 2nd-pointers\n", (int) (nt * nr));
      free ((void *) array);
      return NULL;
    }

  for (i = 1; i < nt; i++)
    {
      array[i] = (void **) ((unsigned char *) array[0] + i * nr * sizeof (void *));
    }

  if ((array[0][0] = (void *) calloc (nt * nr * nc, size)) == NULL)
    {
      printf ("[calloc_3d] failed to alloc. memory (%d X %d X %d of size %d)\n", 
      (int) nt, (int) nr, (int) nc, (int) size);
      free ((void *) array[0]);
      free ((void *) array);
      return NULL;
    }
  for (j = 1; j < nr; j++)
    {
      array[0][j] = (void **) ((unsigned char *) array[0][j - 1] + nc * size);
    }

  for (i = 1; i < nt; i++)
    {
      array[i][0] = (void **) ((unsigned char *) array[i - 1][0] + nr * nc * size);
      for (j = 1; j < nr; j++)
	{
	  array[i][j] = (void **) ((unsigned char *) array[i][j - 1] + nc * size);
	}
    }

  return array;
}


/*----------------------------------------------------------------------------*/
/*! \fn void free_1d_array(void *array)
 *  *  \brief Free memory used by 1D array  */
void
free_1d_array (void *array)
{
  free (array);
}

/*----------------------------------------------------------------------------*/
/*! \fn void free_2d_array(void *array)
 *  *  \brief Free memory used by 2D array  */
void
free_2d_array (void *array)
{
  void **ta = (void **) array;

  free (ta[0]);
  free (array);
}

/*----------------------------------------------------------------------------*/
/*! \fn void free_3d_array(void *array)
 *  *  \brief Free memory used by 3D array  */
void
free_3d_array (void *array)
{
  void ***ta = (void ***) array;

  free (ta[0][0]);
  free (ta[0]);
  free (array);
}
