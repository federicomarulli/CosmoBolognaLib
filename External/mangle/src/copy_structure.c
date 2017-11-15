/*------------------------------------------------------------------------------
© A J S Hamilton 2001
------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------
  Copy structure from struct1 to struct2.

   Input: n = sizeof(struct1).
	  struct1, structure = pointers to structures of same size.
  Output: revised contents of struct2.
*/
void copy_structure(int n, char *struct1, char *struct2)
{
    char *p1, *p2;
    int i;

    p1 = (char *) struct1;
    p2 = (char *) struct2;
    for (i = 0; i < n; i++) *p2++ = *p1++;
}
