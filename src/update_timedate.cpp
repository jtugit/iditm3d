/*************************************************************************
  Function update_time_date_year
 * With advance of the time, universal time sec changes and date, year may also
 * need to be updated

  Input: parms - structure for system parameters
 *       t     - time in sec relative to the start of the simulation
 * Output: parms. change on time, date, & year are stored in parms

  By Jiannan Tu
  1/13/2014
*************************************************************************/
#include "param.h"

void update_timedate(AppCtx *parms)
{
    const  int mons[12]={31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 20, 31};
    int    date, imon, yr;
    double secs;

    imon=parms->mon;
    yr=parms->iyr;
    date=parms->idate+1;
    secs=parms->sec;

    if (secs >= 86400.0) {
        secs=secs-86400.0;

        if (date > mons[imon-1]) {
            switch (imon) {
                case 2:
                    if (yr % 4 == 0) {
                        if (date > 29) {
                            imon += 1;
                            date=1;
                        }
                    }
                    else {
                        imon += 1;
                        date=1;
                    }
                    break;
                case 12:
                    yr += 1;
                    imon=1;
                    date=1;
                    break;
                default:
                    imon += 1;
                    date=1;
            }
        }
    }

    parms->sec=secs;
    parms->idate=date;
    parms->mon=imon;
    parms->iyr=yr;

    return;
}
