#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <libgen.h>

unsigned long dist(char * s1, char * s2, int len);

int main(int argc, char * argv[])
{
    if (argc != 2)
    {
        printf("must provide exactly 1 .data file");
        return -1;
    }

    FILE * f1 = fopen(argv[1], "r");

    char buf[256];
    char * str;

    char ** strings = (char **) malloc(sizeof(char *) * 500000);

    int len = 0;

    int i, ii;

    int istring = 0;
    while (fgets(buf, 256, f1) != NULL)
    {
        if (buf[0] == '#' || buf[0] == ' ' || buf[0] == '\t' || buf[0] == '\n')
            continue;

        str = (char *) malloc(15);

        for (i = 0; buf[i] != '\t' ; i++)
        {
            str[i] = buf[i];
        }
        
        // we save the length of the motifs,
        // so we don't have to call strlen later on
        if (len == 0)
            len = i;

        // null-terminating the string to avoid troubles,
        // but technically we don't need it (we use the above len)
        str[i] = '\0';

        strings[istring] = str;
        istring++;
    }

    if (istring == 0)
    {
        printf("-");
        return 0;
    }
    else if (istring == 1)
    {
        printf("0");
        return 0;
    }

    unsigned long tot = 0;
    unsigned long n = 0;
    for (i = 0; i < istring; i++)
    {
        for (ii = i+1; ii < istring; ii++)
        {
            tot += dist(strings[i], strings[ii], len);

            n++;
        }
    }

    printf("%f", ((double) tot) / (double)(len * n));

    fclose(f1);

    return 0;
}

unsigned long dist(char * s1, char * s2, int len)
{
    unsigned long count = 0;
    
    for (int i = 0; i < len; i++)
    {
        if (s1[i] != s2[i] && s1[i] != 'N' && s2[i] != 'N')
            count++;
    }

    return count;
}