/*************************************************/
/*FileName: test.c                               */
/*Author  : xfeng                                */
/*Mail    : geophydogvon@gmail.com               */
/*Inst    : NJU                                  */
/*Time    : 2017-09-10                           */
/*This is a c programming language!              */
/*************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sacio.h"

int main( int argc, char *argv[] ) {
    char buff[500], sac1[50], sac2[50], cor_name[50], sacbp1[100], sacbp2[100], saccut1[100],
         saccut2[100], sacnorm1[100], sacnorm2[100], sacwhi1[100], sacwhi2[100];
    int year, mon, day, jday, hour, min, sec, cut_npts, npts, count = 1;
    float evt0, f1, f2, f3, f4, lag_time, start0;
    FILE *ff;

    if ( argc != 2 ) {
        fprintf(stderr, "Usage: run_exam file.lst\n");
        exit(1);
    }

    ff = fopen( argv[1], "r" );
    while ( fgets( buff, 500, ff ) ) {
        sscanf(buff, "%s %s %d %d %d %d %d %d %f %d %f %f %f %f %d %s %f", sac1, sac2, &year, &mon, &day,\
            &hour, &min, &sec, &start0, &cut_npts, &f1, &f2, &f3, &f4, &npts, cor_name, &lag_time );

        strcpy(saccut1, sac1);   strcpy(saccut2, sac2);   strcpy(sacbp1, sac1);  strcpy(sacbp2, sac2);
        strcat(saccut1, ".cut"); strcat(saccut2, ".cut"); strcat(sacbp1, ".bp"); strcat(sacbp2, ".bp");

        strcpy(sacnorm1, sac1);    strcpy(sacnorm2, sac2);    strcpy(sacwhi1, sac1);   strcpy(sacwhi2, sac2);
        strcat(sacnorm1, ".norm"); strcat(sacnorm2, ".norm"); strcat(sacwhi1, ".whi"); strcat(sacwhi2, ".whi");

        if ( count % 50 == 0 ) printf("%d\n", count);
        jday = julian(year, mon, day);
        evt0 = abs_time( year, jday, hour, min, sec, 0. );
        cut_sac( sac1, saccut1, evt0, start0, cut_npts );
        bp( saccut1, sacbp1, f1, f2, f3, f4, 10);
        normal( sacbp1, sacnorm1, npts );
        spe_whi ( sacnorm1, sacwhi1, 20, f1, f2, f3, f4 );

        cut_sac( sac2, saccut2, evt0, start0, cut_npts );
        bp( saccut2, sacbp2, f1, f2, f3, f4, 10);
        normal( sacbp2, sacnorm2, npts );
        spe_whi ( sacnorm2, sacwhi2, 20, f1, f2, f3, f4 );

        cor_in_freq( sacwhi1, sacwhi2, lag_time, cor_name );
        count += 1;
    }
    fclose(ff);
    system("mkdir COR"); system("mv COR*.SAC COR/"); system("rm *.cut *.bp *.norm *.whi");
    system("mv COR ../");

    return 0;
}
