#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>

int main(int argc, char * argv[])
{
    if (argc < 4)
        return -1;

    FILE * f1 = fopen(argv[1], "r");
    FILE * f2 = fopen(argv[2], "r");
    FILE * f3 = fopen(argv[3], "w");

    char str1[15];
    char str2[15];

    char ** f1_strs = (char **) malloc(sizeof(char *) * 500000);
    char ** f2_strs = (char **) malloc(sizeof(char *) * 500000);

    int len_f1 = 0;
    int len_f2 = 0;

    int i,ii,iii;

    printf("Loading file: %s\n", argv[1]);
    while (fgets(str1, 12, f1) != NULL)
    {
        f1_strs[len_f1] = (char *) malloc(sizeof(char) * 15);
        str1[10] = '\0';
        strncpy(f1_strs[len_f1], str1, 11);
        len_f1++;
    }

    printf("Loading file: %s\n", argv[2]);
    while (fgets(str2, 12, f2) != NULL)
    {
        f2_strs[len_f2] = (char *) malloc(sizeof(char) * 15);
        str2[10] = '\0';
        strncpy(f2_strs[len_f2], str2, 11);
        len_f2++;
    }

    for (i = 0; i < len_f1; i++)
    {
        for (ii = 0; ii < len_f2; ii++)
        {
            for (iii = 0; iii < 10; iii++)
            {
                if (f1_strs[i][iii] != f2_strs[ii][iii] && f1_strs[i][iii] != 'N' && f2_strs[ii][iii] != 'N')
                    break;
            }
            if (iii == 10)
            {
                f1_strs[i][10] = '\n';
                f2_strs[ii][10] = '\n';
                fputs(f1_strs[i], f3);
                fputs(f2_strs[ii], f3);
            }
        }
    }

    return 0;
}
