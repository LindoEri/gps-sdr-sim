using System;
using System.Text;
using System.Threading;
using static System.Math;
using System.IO;
using System.Diagnostics;

namespace GPS模拟
{
    //    MIT License

    //Copyright(c) 2015 Takuji Ebinuma

    //Permission is hereby granted, free of charge, to any person obtaining a copy
    //of this software and associated documentation files (the "Software"), to deal
    //in the Software without restriction, including without limitation the rights
    //to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    //copies of the Software, and to permit persons to whom the Software is
    //furnished to do so, subject to the following conditions:

    //The above copyright notice and this permission notice shall be included in all
    //copies or substantial portions of the Software.

    //THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    //IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    //FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    //AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    //LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    //OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    //SOFTWARE.


    static class GPS
    {

        const int MAX_CHAR = 100;
        const int MAX_SAT = 32;
        const int MAX_CHAN = 16;
        const int USER_MOTION_SIZE = 3000;
        const int STATIC_MAX_DURATION = 86400;
        const int N_SBF = 5;
        const int N_DWRD_SBF = 10;
        const int N_DWRD = (N_SBF + 1) * N_DWRD_SBF;
        const int CA_SEQ_LEN = 1023;
        const double SECONDS_IN_WEEK = 604800.0;
        const double SECONDS_IN_HALF_WEEK = 302400.0;
        const double SECONDS_IN_DAY = 86400.0;
        const double SECONDS_IN_HOUR = 3600.0;
        const double SECONDS_IN_MINUTE = 60.0;
        const double POW2_M5 = 0.03125;
        const double POW2_M19 = 1.907348632812500e-6;
        const double POW2_M29 = 1.862645149230957e-9;
        const double POW2_M31 = 4.656612873077393e-10;
        const double POW2_M33 = 1.164153218269348e-10;
        const double POW2_M43 = 1.136868377216160e-13;
        const double POW2_M55 = 2.775557561562891e-17;
        const double POW2_M50 = 8.881784197001252e-016;
        const double POW2_M30 = 9.313225746154785e-010;
        const double POW2_M27 = 7.450580596923828e-009;
        const double POW2_M24 = 5.960464477539063e-008;
        const double GM_EARTH = 3.986005e14;
        const double OMEGA_EARTH = 7.2921151467e-5;
        const double PI = 3.1415926535898;
        const double WGS84_RADIUS = 6378137.0;
        const double WGS84_ECCENTRICITY = 0.0818191908426;
        const double R2D = 57.2957795131;
        const double SPEED_OF_LIGHT = 2.99792458e8;
        const double LAMBDA_L1 = 0.190293672798365;
        const double CARR_FREQ = 1575.42e6;
        const double CODE_FREQ = 1.023e6;
        const double CARR_TO_CODE = 1.0 / 1540.0;
        const int SC01 = 1;
        const int SC08 = 8;
        const int SC16 = 16;
        const int EPHEM_ARRAY_SIZE = 13;

        struct gpstime_t
        {
            public int week;   /*!< GPS week number (since January 1980) */
            public double sec;     /*!< second inside the GPS \a week */
        };
        struct datetime_t
        {
            public int y;      /*!< Calendar year */
            public int m;      /*!< Calendar month */
            public int d;      /*!< Calendar day */
            public int hh;     /*!< Calendar hour */
            public int mm;     /*!< Calendar minutes */
            public double sec; /*!< Calendar seconds */
        };
        struct ephem_t
        {
            public int vflg;   /*!< Valid Flag */
            public datetime_t t;
            public gpstime_t toc;  /*!< Time of Clock */
            public gpstime_t toe;  /*!< Time of Ephemeris */
            public int iodc;   /*!< Issue of Data, Clock */
            public int iode;   /*!< Isuse of Data, Ephemeris */
            public double deltan;  /*!< Delta-N (radians/sec) */
            public double cuc; /*!< Cuc (radians) */
            public double cus; /*!< Cus (radians) */
            public double cic; /*!< Correction to inclination cos (radians) */
            public double cis; /*!< Correction to inclination sin (radians) */
            public double crc; /*!< Correction to radius cos (meters) */
            public double crs; /*!< Correction to radius sin (meters) */
            public double ecc; /*!< e Eccentricity */
            public double sqrta;   /*!< sqrt(A) (sqrt(m)) */
            public double m0;  /*!< Mean anamoly (radians) */
            public double omg0;    /*!< Longitude of the ascending node (radians) */
            public double inc0;    /*!< Inclination (radians) */
            public double aop;
            public double omgdot;  /*!< Omega dot (radians/s) */
            public double idot;    /*!< IDOT (radians/s) */
            public double af0; /*!< Clock offset (seconds) */
            public double af1; /*!< rate (sec/sec) */
            public double af2; /*!< acceleration (sec/sec^2) */
            public double tgd; /*!< Group delay L2 bias */
            public int svhlth;
            public int codeL2;
            // Working variables follow
            public double n;   /*!< Mean motion (Average angular velocity) */
            public double sq1e2;   /*!< sqrt(1-e^2) */
            public double A;   /*!< Semi-major axis */
            public double omgkdot; /*!< OmegaDot-OmegaEdot */
        };
        struct ionoutc_t
        {

            public int enable;
            public int vflg;
            public double alpha0, alpha1, alpha2, alpha3;
            public double beta0, beta1, beta2, beta3;
            public double A0, A1;
            public int dtls, tot, wnt;
            public int dtlsf, dn, wnlsf;
        };
        struct range_t
        {

            public gpstime_t g;
            public double range; // pseudorange
            public double rate;
            public double d; // geometric distance
            public double[] azel;
            public double iono_delay;
        };
        struct channel_t
        {

            public int prn;    /*< PRN Number */
            public int[] ca; /*< C/A Sequence */
            public double f_carr;  /*< Carrier frequency */
            public double f_code;  /*< Code frequency */
            public double carr_phase;

            public double code_phase; /*< Code phase */
            public gpstime_t g0;   /*!< GPS time at start */
            public ulong[][] sbf; /*!< current subframe */
            public ulong[] dwrd; /*!< Data words of sub-frame */
            public int iword;  /*!< initial word */
            public int ibit;   /*!< initial bit */
            public int icode;  /*!< initial code */
            public int dataBit;    /*!< current data bit */
            public int codeCA; /*!< current C/A code */
            public double[] azel;
            public range_t rho0;
            public int carr_phasestep;  /*< Carrier phasestep */
        };

        static int[] sinTable512 = {
    2,   5,   8,  11,  14,  17,  20,  23,  26,  29,  32,  35,  38,  41,  44,  47,
    50,  53,  56,  59,  62,  65,  68,  71,  74,  77,  80,  83,  86,  89,  91,  94,
    97, 100, 103, 105, 108, 111, 114, 116, 119, 122, 125, 127, 130, 132, 135, 138,
    140, 143, 145, 148, 150, 153, 155, 157, 160, 162, 164, 167, 169, 171, 173, 176,
    178, 180, 182, 184, 186, 188, 190, 192, 194, 196, 198, 200, 202, 204, 205, 207,
    209, 210, 212, 214, 215, 217, 218, 220, 221, 223, 224, 225, 227, 228, 229, 230,
    232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 241, 242, 243, 244, 244, 245,
    245, 246, 247, 247, 248, 248, 248, 249, 249, 249, 249, 250, 250, 250, 250, 250,
    250, 250, 250, 250, 250, 249, 249, 249, 249, 248, 248, 248, 247, 247, 246, 245,
    245, 244, 244, 243, 242, 241, 241, 240, 239, 238, 237, 236, 235, 234, 233, 232,
    230, 229, 228, 227, 225, 224, 223, 221, 220, 218, 217, 215, 214, 212, 210, 209,
    207, 205, 204, 202, 200, 198, 196, 194, 192, 190, 188, 186, 184, 182, 180, 178,
    176, 173, 171, 169, 167, 164, 162, 160, 157, 155, 153, 150, 148, 145, 143, 140,
    138, 135, 132, 130, 127, 125, 122, 119, 116, 114, 111, 108, 105, 103, 100,  97,
    94,  91,  89,  86,  83,  80,  77,  74,  71,  68,  65,  62,  59,  56,  53,  50,
    47,  44,  41,  38,  35,  32,  29,  26,  23,  20,  17,  14,  11,   8,   5,   2,
    -2,  -5,  -8, -11, -14, -17, -20, -23, -26, -29, -32, -35, -38, -41, -44, -47,
    -50, -53, -56, -59, -62, -65, -68, -71, -74, -77, -80, -83, -86, -89, -91, -94,
    -97,-100,-103,-105,-108,-111,-114,-116,-119,-122,-125,-127,-130,-132,-135,-138,
    -140,-143,-145,-148,-150,-153,-155,-157,-160,-162,-164,-167,-169,-171,-173,-176,
    -178,-180,-182,-184,-186,-188,-190,-192,-194,-196,-198,-200,-202,-204,-205,-207,
    -209,-210,-212,-214,-215,-217,-218,-220,-221,-223,-224,-225,-227,-228,-229,-230,
    -232,-233,-234,-235,-236,-237,-238,-239,-240,-241,-241,-242,-243,-244,-244,-245,
    -245,-246,-247,-247,-248,-248,-248,-249,-249,-249,-249,-250,-250,-250,-250,-250,
    -250,-250,-250,-250,-250,-249,-249,-249,-249,-248,-248,-248,-247,-247,-246,-245,
    -245,-244,-244,-243,-242,-241,-241,-240,-239,-238,-237,-236,-235,-234,-233,-232,
    -230,-229,-228,-227,-225,-224,-223,-221,-220,-218,-217,-215,-214,-212,-210,-209,
    -207,-205,-204,-202,-200,-198,-196,-194,-192,-190,-188,-186,-184,-182,-180,-178,
    -176,-173,-171,-169,-167,-164,-162,-160,-157,-155,-153,-150,-148,-145,-143,-140,
    -138,-135,-132,-130,-127,-125,-122,-119,-116,-114,-111,-108,-105,-103,-100, -97,
    -94, -91, -89, -86, -83, -80, -77, -74, -71, -68, -65, -62, -59, -56, -53, -50,
    -47, -44, -41, -38, -35, -32, -29, -26, -23, -20, -17, -14, -11,  -8,  -5,  -2};
        static int[] cosTable512 = {
    250, 250, 250, 250, 250, 249, 249, 249, 249, 248, 248, 248, 247, 247, 246, 245,
    245, 244, 244, 243, 242, 241, 241, 240, 239, 238, 237, 236, 235, 234, 233, 232,
    230, 229, 228, 227, 225, 224, 223, 221, 220, 218, 217, 215, 214, 212, 210, 209,
    207, 205, 204, 202, 200, 198, 196, 194, 192, 190, 188, 186, 184, 182, 180, 178,
    176, 173, 171, 169, 167, 164, 162, 160, 157, 155, 153, 150, 148, 145, 143, 140,
    138, 135, 132, 130, 127, 125, 122, 119, 116, 114, 111, 108, 105, 103, 100,  97,
    94,  91,  89,  86,  83,  80,  77,  74,  71,  68,  65,  62,  59,  56,  53,  50,
    47,  44,  41,  38,  35,  32,  29,  26,  23,  20,  17,  14,  11,   8,   5,   2,
    -2,  -5,  -8, -11, -14, -17, -20, -23, -26, -29, -32, -35, -38, -41, -44, -47,
    -50, -53, -56, -59, -62, -65, -68, -71, -74, -77, -80, -83, -86, -89, -91, -94,
    -97,-100,-103,-105,-108,-111,-114,-116,-119,-122,-125,-127,-130,-132,-135,-138,
    -140,-143,-145,-148,-150,-153,-155,-157,-160,-162,-164,-167,-169,-171,-173,-176,
    -178,-180,-182,-184,-186,-188,-190,-192,-194,-196,-198,-200,-202,-204,-205,-207,
    -209,-210,-212,-214,-215,-217,-218,-220,-221,-223,-224,-225,-227,-228,-229,-230,
    -232,-233,-234,-235,-236,-237,-238,-239,-240,-241,-241,-242,-243,-244,-244,-245,
    -245,-246,-247,-247,-248,-248,-248,-249,-249,-249,-249,-250,-250,-250,-250,-250,
    -250,-250,-250,-250,-250,-249,-249,-249,-249,-248,-248,-248,-247,-247,-246,-245,
    -245,-244,-244,-243,-242,-241,-241,-240,-239,-238,-237,-236,-235,-234,-233,-232,
    -230,-229,-228,-227,-225,-224,-223,-221,-220,-218,-217,-215,-214,-212,-210,-209,
    -207,-205,-204,-202,-200,-198,-196,-194,-192,-190,-188,-186,-184,-182,-180,-178,
    -176,-173,-171,-169,-167,-164,-162,-160,-157,-155,-153,-150,-148,-145,-143,-140,
    -138,-135,-132,-130,-127,-125,-122,-119,-116,-114,-111,-108,-105,-103,-100, -97,
    -94, -91, -89, -86, -83, -80, -77, -74, -71, -68, -65, -62, -59, -56, -53, -50,
    -47, -44, -41, -38, -35, -32, -29, -26, -23, -20, -17, -14, -11,  -8,  -5,  -2,
    2,   5,   8,  11,  14,  17,  20,  23,  26,  29,  32,  35,  38,  41,  44,  47,
    50,  53,  56,  59,  62,  65,  68,  71,  74,  77,  80,  83,  86,  89,  91,  94,
    97, 100, 103, 105, 108, 111, 114, 116, 119, 122, 125, 127, 130, 132, 135, 138,
    140, 143, 145, 148, 150, 153, 155, 157, 160, 162, 164, 167, 169, 171, 173, 176,
    178, 180, 182, 184, 186, 188, 190, 192, 194, 196, 198, 200, 202, 204, 205, 207,
    209, 210, 212, 214, 215, 217, 218, 220, 221, 223, 224, 225, 227, 228, 229, 230,
    232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 241, 242, 243, 244, 244, 245,
    245, 246, 247, 247, 248, 248, 248, 249, 249, 249, 249, 250, 250, 250, 250, 250};
        static double[] ant_pat_db = {
    0.00,  0.00,  0.22,  0.44,  0.67,  1.11,  1.56,  2.00,  2.44,  2.89,  3.56,  4.22,
    4.89,  5.56,  6.22,  6.89,  7.56,  8.22,  8.89,  9.78, 10.67, 11.56, 12.44, 13.33,
    14.44, 15.56, 16.67, 17.78, 18.89, 20.00, 21.33, 22.67, 24.00, 25.56, 27.33, 29.33,
    31.56};
        static int[] allocatedSat = new int[32];

        static void subVect(double[] y, double[] x1, double[] x2)
        {
            y[0] = x1[0] - x2[0];
            y[1] = x1[1] - x2[1];
            y[2] = x1[2] - x2[2];
            return;
        }
        static double normVect(double[] x)
        {
            return (Sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]));
        }
        static double dotProd(double[] x1, double[] x2)
        {
            return (x1[0] * x2[0] + x1[1] * x2[1] + x1[2] * x2[2]);
        }
        static void codegen(int[] ca, int prn)
        {
            int[] delay = {
        5,   6,   7,   8,  17,  18, 139, 140, 141, 251,
        252, 254, 255, 256, 257, 258, 469, 470, 471, 472,
        473, 474, 509, 512, 513, 514, 515, 516, 859, 860,
        861, 862 };

            int[] g1 = new int[CA_SEQ_LEN], g2 = new int[CA_SEQ_LEN];
            int[] r1 = new int[N_DWRD_SBF], r2 = new int[N_DWRD_SBF];
            int c1, c2;
            int i, j;

            if (prn < 1 || prn > 32)
                return;

            for (i = 0; i < N_DWRD_SBF; i++)
                r1[i] = r2[i] = -1;

            for (i = 0; i < CA_SEQ_LEN; i++)
            {
                g1[i] = r1[9];
                g2[i] = r2[9];
                c1 = r1[2] * r1[9];
                c2 = r2[1] * r2[2] * r2[5] * r2[7] * r2[8] * r2[9];

                for (j = 9; j > 0; j--)
                {
                    r1[j] = r1[j - 1];
                    r2[j] = r2[j - 1];
                }
                r1[0] = c1;
                r2[0] = c2;
            }

            for (i = 0, j = CA_SEQ_LEN - delay[prn - 1]; i < CA_SEQ_LEN; i++, j++)
                ca[i] = (1 - g1[i] * g2[j % CA_SEQ_LEN]) / 2;

            return;
        }
        static void date2gps(datetime_t t, ref gpstime_t g)
        {
            int[] doy = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 };
            int ye;
            int de;
            int lpdays;

            ye = t.y - 1980;

            // Compute the number of leap days since Jan 5/Jan 6, 1980.
            lpdays = ye / 4 + 1;
            if ((ye % 4) == 0 && t.m <= 2)
                lpdays--;

            // Compute the number of days elapsed since Jan 5/Jan 6, 1980.
            de = ye * 365 + doy[t.m - 1] + t.d + lpdays - 6;

            // Convert time to GPS weeks and seconds.
            g.week = de / 7;
            g.sec = (double)(de % 7) * 86400 + t.hh * 3600 + t.mm * 60 + t.sec;

            return;
        }
        static void gps2date(gpstime_t g, ref datetime_t t)
        {
            // Convert Julian day number to calendar date
            int c = (int)(7 * g.week + Floor(g.sec / 86400.0) + 2444245.0) + 1537;
            int d = (int)((c - 122.1) / 365.25);
            int e = 365 * d + d / 4;
            int f = (int)((c - e) / 30.6001);

            t.d = c - e - (int)(30.6001 * f);
            t.m = f - 1 - 12 * (f / 14);
            t.y = d - 4715 - ((7 + t.m) / 10);

            t.hh = ((int)(g.sec / 3600.0)) % 24;
            t.mm = ((int)(g.sec / 60.0)) % 60;
            t.sec = g.sec - 60.0 * Floor(g.sec / 60.0);

            return;
        }
        static void xyz2llh(double[] xyz, double[] llh)
        {

            double a, eps, e, e2;
            double x, y, z;
            double rho2, dz, zdz, nh, slat, n, dz_new;

            a = WGS84_RADIUS;
            e = WGS84_ECCENTRICITY;

            eps = 1.0e-3;
            e2 = e * e;

            if (normVect(xyz) < eps)
            {
                // Invalid ECEF vector
                llh[0] = 0.0;
                llh[1] = 0.0;
                llh[2] = -a;

                return;
            }

            x = xyz[0];
            y = xyz[1];
            z = xyz[2];

            rho2 = x * x + y * y;
            dz = e2 * z;

            while (true)
            {
                zdz = z + dz;
                nh = Sqrt(rho2 + zdz * zdz);
                slat = zdz / nh;
                n = a / Sqrt(1.0 - e2 * slat * slat);
                dz_new = n * e2 * slat;

                if (Abs(dz - dz_new) < eps)
                    break;

                dz = dz_new;
            }

            llh[0] = Atan2(zdz, Sqrt(rho2));
            llh[1] = Atan2(y, x);
            llh[2] = nh - n;

            return;
        }
        static void llh2xyz(double[] llh, double[] xyz)
        {
            double n;
            double a;
            double e;
            double e2;
            double clat;
            double slat;
            double clon;
            double slon;
            double d, nph;
            double tmp;

            a = WGS84_RADIUS;
            e = WGS84_ECCENTRICITY;
            e2 = e * e;

            clat = Cos(llh[0]);
            slat = Sin(llh[0]);
            clon = Cos(llh[1]);
            slon = Sin(llh[1]);
            d = e * slat;

            n = a / Sqrt(1.0 - d * d);
            nph = n + llh[2];

            tmp = nph * clat;

            double rad = 0, phi = 0, theta = 0;
            if (inacc != 0)
            {
                Random r = new Random();
                rad = (r.NextDouble() - 0.5) * 2 * inacc;
                phi = (r.NextDouble() - 0.5) * 2 * 90;
                theta = r.NextDouble() * 360;
            }


            xyz[0] = tmp * clon + rad * Sin(phi) * Sin(theta);
            xyz[1] = tmp * slon + rad * Sin(phi) * Cos(theta);
            xyz[2] = ((1.0 - e2) * n + llh[2]) * slat + rad * Cos(phi);

            return;
        }
        static void ltcmat(double[] llh, double[][] t)
        {
            double slat, clat;
            double slon, clon;

            slat = Sin(llh[0]);
            clat = Cos(llh[0]);
            slon = Sin(llh[1]);
            clon = Cos(llh[1]);

            t[0][0] = -slat * clon;
            t[0][1] = -slat * slon;
            t[0][2] = clat;
            t[1][0] = -slon;
            t[1][1] = clon;
            t[1][2] = 0.0;
            t[2][0] = clat * clon;
            t[2][1] = clat * slon;
            t[2][2] = slat;

            return;
        }
        static void ecef2neu(double[] xyz, double[][] t, double[] neu)
        {
            neu[0] = t[0][0] * xyz[0] + t[0][1] * xyz[1] + t[0][2] * xyz[2];
            neu[1] = t[1][0] * xyz[0] + t[1][1] * xyz[1] + t[1][2] * xyz[2];
            neu[2] = t[2][0] * xyz[0] + t[2][1] * xyz[1] + t[2][2] * xyz[2];

            return;
        }
        static void neu2azel(double[] azel, double[] neu)
        {
            double ne;

            azel[0] = Atan2(neu[1], neu[0]);
            if (azel[0] < 0.0)
                azel[0] += (2.0 * PI);

            ne = Sqrt(neu[0] * neu[0] + neu[1] * neu[1]);
            azel[1] = Atan2(neu[2], ne);

            return;
        }
        static void satpos(ephem_t eph, gpstime_t g, double[] pos, double[] vel, double[] clk)
        {

            double tk;
            double mk;
            double ek;
            double ekold;
            double ekdot;
            double cek, sek;
            double pk;
            double pkdot;
            double c2pk, s2pk;
            double uk;
            double ukdot;
            double cuk, suk;
            double ok;
            double sok, cok;
            double ik;
            double ikdot;
            double sik, cik;
            double rk;
            double rkdot;
            double xpk, ypk;
            double xpkdot, ypkdot;

            double relativistic, OneMinusecosE, tmp;

            tk = g.sec - eph.toe.sec;

            if (tk > SECONDS_IN_HALF_WEEK)
                tk -= SECONDS_IN_WEEK;
            else if (tk < -SECONDS_IN_HALF_WEEK)
                tk += SECONDS_IN_WEEK;

            mk = eph.m0 + eph.n * tk;
            ek = mk;
            ekold = ek + 1.0;

            OneMinusecosE = 0; // Suppress the uninitialized warning.
            while (Abs(ek - ekold) > 1.0E-14)
            {
                ekold = ek;
                OneMinusecosE = 1.0 - eph.ecc * Cos(ekold);
                ek = ek + (mk - ekold + eph.ecc * Sin(ekold)) / OneMinusecosE;
            }

            sek = Sin(ek);
            cek = Cos(ek);

            ekdot = eph.n / OneMinusecosE;

            relativistic = -4.442807633E-10 * eph.ecc * eph.sqrta * sek;

            pk = Atan2(eph.sq1e2 * sek, cek - eph.ecc) + eph.aop;
            pkdot = eph.sq1e2 * ekdot / OneMinusecosE;

            s2pk = Sin(2.0 * pk);
            c2pk = Cos(2.0 * pk);

            uk = pk + eph.cus * s2pk + eph.cuc * c2pk;
            suk = Sin(uk);
            cuk = Cos(uk);
            ukdot = pkdot * (1.0 + 2.0 * (eph.cus * c2pk - eph.cuc * s2pk));

            rk = eph.A * OneMinusecosE + eph.crc * c2pk + eph.crs * s2pk;
            rkdot = eph.A * eph.ecc * sek * ekdot + 2.0 * pkdot * (eph.crs * c2pk - eph.crc * s2pk);

            ik = eph.inc0 + eph.idot * tk + eph.cic * c2pk + eph.cis * s2pk;
            sik = Sin(ik);
            cik = Cos(ik);
            ikdot = eph.idot + 2.0 * pkdot * (eph.cis * c2pk - eph.cic * s2pk);

            xpk = rk * cuk;
            ypk = rk * suk;
            xpkdot = rkdot * cuk - ypk * ukdot;
            ypkdot = rkdot * suk + xpk * ukdot;

            ok = eph.omg0 + tk * eph.omgkdot - OMEGA_EARTH * eph.toe.sec;
            sok = Sin(ok);
            cok = Cos(ok);

            pos[0] = xpk * cok - ypk * cik * sok;
            pos[1] = xpk * sok + ypk * cik * cok;
            pos[2] = ypk * sik;

            tmp = ypkdot * cik - ypk * sik * ikdot;

            vel[0] = -eph.omgkdot * pos[1] + xpkdot * cok - tmp * sok;
            vel[1] = eph.omgkdot * pos[0] + xpkdot * sok + tmp * cok;
            vel[2] = ypk * cik * ikdot + ypkdot * sik;

            // Satellite clock correction
            tk = g.sec - eph.toc.sec;

            if (tk > SECONDS_IN_HALF_WEEK)
                tk -= SECONDS_IN_WEEK;
            else if (tk < -SECONDS_IN_HALF_WEEK)
                tk += SECONDS_IN_WEEK;

            clk[0] = eph.af0 + tk * (eph.af1 + tk * eph.af2) + relativistic - eph.tgd;
            clk[1] = eph.af1 + 2.0 * tk * eph.af2;

            return;
        }
        static void eph2sbf(ephem_t eph, ionoutc_t ionoutc, ulong[][] sbf)
        {
            ulong wn;
            ulong toe;
            ulong toc;
            ulong iode;
            ulong iodc;
            long deltan;
            long cuc;
            long cus;
            long cic;
            long cis;
            long crc;
            long crs;
            ulong ecc;
            ulong sqrta;
            long m0;
            long omg0;
            long inc0;
            long aop;
            long omgdot;
            long idot;
            long af0;
            long af1;
            long af2;
            long tgd;
            int svhlth;
            int codeL2;

            ulong ura = 0UL;
            ulong dataId = 1UL;
            ulong sbf4_page25_svId = 63UL;
            ulong sbf5_page25_svId = 51UL;

            ulong wna;
            ulong toa;

            long alpha0, alpha1, alpha2, alpha3;
            long beta0, beta1, beta2, beta3;
            long A0, A1;
            long dtls, dtlsf;
            ulong tot, wnt, wnlsf, dn;
            ulong sbf4_page18_svId = 56UL;

            // FIXED: This has to be the "transmission" week number, not for the ephemeris reference time
            //wn = (ulong)(eph.toe.week%1024);
            wn = 0UL;
            toe = (ulong)(eph.toe.sec / 16.0);
            toc = (ulong)(eph.toc.sec / 16.0);
            iode = (ulong)(eph.iode);
            iodc = (ulong)(eph.iodc);
            deltan = (long)(eph.deltan / POW2_M43 / PI);
            cuc = (long)(eph.cuc / POW2_M29);
            cus = (long)(eph.cus / POW2_M29);
            cic = (long)(eph.cic / POW2_M29);
            cis = (long)(eph.cis / POW2_M29);
            crc = (long)(eph.crc / POW2_M5);
            crs = (long)(eph.crs / POW2_M5);
            ecc = (ulong)(eph.ecc / POW2_M33);
            sqrta = (ulong)(eph.sqrta / POW2_M19);
            m0 = (long)(eph.m0 / POW2_M31 / PI);
            omg0 = (long)(eph.omg0 / POW2_M31 / PI);
            inc0 = (long)(eph.inc0 / POW2_M31 / PI);
            aop = (long)(eph.aop / POW2_M31 / PI);
            omgdot = (long)(eph.omgdot / POW2_M43 / PI);
            idot = (long)(eph.idot / POW2_M43 / PI);
            af0 = (long)(eph.af0 / POW2_M31);
            af1 = (long)(eph.af1 / POW2_M43);
            af2 = (long)(eph.af2 / POW2_M55);
            tgd = (long)(eph.tgd / POW2_M31);
            svhlth = (int)(eph.svhlth);
            codeL2 = (int)(eph.codeL2);

            wna = (ulong)(eph.toe.week % 256);
            toa = (ulong)(eph.toe.sec / 4096.0);

            alpha0 = (long)Round(ionoutc.alpha0 / POW2_M30);
            alpha1 = (long)Round(ionoutc.alpha1 / POW2_M27);
            alpha2 = (long)Round(ionoutc.alpha2 / POW2_M24);
            alpha3 = (long)Round(ionoutc.alpha3 / POW2_M24);
            beta0 = (long)Round(ionoutc.beta0 / 2048.0);
            beta1 = (long)Round(ionoutc.beta1 / 16384.0);
            beta2 = (long)Round(ionoutc.beta2 / 65536.0);
            beta3 = (long)Round(ionoutc.beta3 / 65536.0);
            A0 = (long)Round(ionoutc.A0 / POW2_M30);
            A1 = (long)Round(ionoutc.A1 / POW2_M50);
            dtls = (long)(ionoutc.dtls);
            tot = (ulong)(ionoutc.tot / 4096);
            wnt = (ulong)(ionoutc.wnt % 256);
            // TO DO: Specify scheduled leap seconds in command options
            // 2016/12/31 (Sat) . WNlsf = 1929, DN = 7 (http://navigationservices.agi.com/GNSSWeb/)
            // Days are counted from 1 to 7 (Sunday is 1).
            wnlsf = 1929 % 256;
            dn = 7;
            dtlsf = 18;

            // Subframe 1
            sbf[0][0] = 0x8B0000UL << 6;
            sbf[0][1] = 0x1UL << 8;
            sbf[0][2] = ((wn & 0x3FFUL) << 20) | (((ulong)codeL2 & 0x3UL) << 18) | ((ura & 0xFUL) << 14) | (((ulong)svhlth & 0x3FUL) << 8) | (((iodc >> 8) & 0x3UL) << 6);
            sbf[0][3] = 0UL;
            sbf[0][4] = 0UL;
            sbf[0][5] = 0UL;
            sbf[0][6] = ((ulong)tgd & 0xFFUL) << 6;
            sbf[0][7] = ((iodc & 0xFFUL) << 22) | ((toc & 0xFFFFUL) << 6);
            sbf[0][8] = (((ulong)af2 & 0xFFUL) << 22) | (((ulong)af1 & 0xFFFFUL) << 6);
            sbf[0][9] = ((ulong)af0 & 0x3FFFFFUL) << 8;

            // Subframe 2
            sbf[1][0] = 0x8B0000UL << 6;
            sbf[1][1] = 0x2UL << 8;
            sbf[1][2] = ((iode & 0xFFUL) << 22) | (((ulong)crs & 0xFFFFUL) << 6);
            sbf[1][3] = (((ulong)deltan & 0xFFFFUL) << 14) | ((((ulong)m0 >> 24) & 0xFFUL) << 6);
            sbf[1][4] = ((ulong)m0 & 0xFFFFFFUL) << 6;
            sbf[1][5] = (((ulong)cuc & 0xFFFFUL) << 14) | (((ecc >> 24) & 0xFFUL) << 6);
            sbf[1][6] = (ecc & 0xFFFFFFUL) << 6;
            sbf[1][7] = (((ulong)cus & 0xFFFFUL) << 14) | (((sqrta >> 24) & 0xFFUL) << 6);
            sbf[1][8] = (sqrta & 0xFFFFFFUL) << 6;
            sbf[1][9] = (toe & 0xFFFFUL) << 14;

            // Subframe 3
            sbf[2][0] = 0x8B0000UL << 6;
            sbf[2][1] = 0x3UL << 8;
            sbf[2][2] = (((ulong)cic & 0xFFFFUL) << 14) | ((((ulong)omg0 >> 24) & 0xFFUL) << 6);
            sbf[2][3] = ((ulong)omg0 & 0xFFFFFFUL) << 6;
            sbf[2][4] = (((ulong)cis & 0xFFFFUL) << 14) | ((((ulong)inc0 >> 24) & 0xFFUL) << 6);
            sbf[2][5] = ((ulong)inc0 & 0xFFFFFFUL) << 6;
            sbf[2][6] = (((ulong)crc & 0xFFFFUL) << 14) | ((((ulong)aop >> 24) & 0xFFUL) << 6);
            sbf[2][7] = ((ulong)aop & 0xFFFFFFUL) << 6;
            sbf[2][8] = ((ulong)omgdot & 0xFFFFFFUL) << 6;
            sbf[2][9] = ((iode & 0xFFUL) << 22) | (((ulong)idot & 0x3FFFUL) << 8);

            if (ionoutc.vflg != 0)
            {
                // Subframe 4, page 18
                sbf[3][0] = 0x8B0000UL << 6;
                sbf[3][1] = 0x4UL << 8;
                sbf[3][2] = (dataId << 28) | (sbf4_page18_svId << 22) | (((ulong)alpha0 & 0xFFUL) << 14) | (((ulong)alpha1 & 0xFFUL) << 6);
                sbf[3][3] = (((ulong)alpha2 & 0xFFUL) << 22) | (((ulong)alpha3 & 0xFFUL) << 14) | (((ulong)beta0 & 0xFFUL) << 6);
                sbf[3][4] = (((ulong)beta1 & 0xFFUL) << 22) | (((ulong)beta2 & 0xFFUL) << 14) | (((ulong)beta3 & 0xFFUL) << 6);
                sbf[3][5] = ((ulong)A1 & 0xFFFFFFUL) << 6;
                sbf[3][6] = (((ulong)A0 >> 8) & 0xFFFFFFUL) << 6;
                sbf[3][7] = (((ulong)A0 & 0xFFUL) << 22) | ((tot & 0xFFUL) << 14) | ((wnt & 0xFFUL) << 6);
                sbf[3][8] = (((ulong)dtls & 0xFFUL) << 22) | ((wnlsf & 0xFFUL) << 14) | ((dn & 0xFFUL) << 6);
                sbf[3][9] = ((ulong)dtlsf & 0xFFUL) << 22;

            }
            else
            {
                // Subframe 4, page 25
                sbf[3][0] = 0x8B0000UL << 6;
                sbf[3][1] = 0x4UL << 8;
                sbf[3][2] = (dataId << 28) | (sbf4_page25_svId << 22);
                sbf[3][3] = 0UL;
                sbf[3][4] = 0UL;
                sbf[3][5] = 0UL;
                sbf[3][6] = 0UL;
                sbf[3][7] = 0UL;
                sbf[3][8] = 0UL;
                sbf[3][9] = 0UL;
            }

            // Subframe 5, page 25
            sbf[4][0] = 0x8B0000UL << 6;
            sbf[4][1] = 0x5UL << 8;
            sbf[4][2] = (dataId << 28) | (sbf5_page25_svId << 22) | ((toa & 0xFFUL) << 14) | ((wna & 0xFFUL) << 6);
            sbf[4][3] = 0UL;
            sbf[4][4] = 0UL;
            sbf[4][5] = 0UL;
            sbf[4][6] = 0UL;
            sbf[4][7] = 0UL;
            sbf[4][8] = 0UL;
            sbf[4][9] = 0UL;

            return;
        }
        static ulong countBits(ulong v)
        {
            ulong c;
            int[] S = { 1, 2, 4, 8, 16 };
            ulong[] B = { 0x55555555, 0x33333333, 0x0F0F0F0F, 0x00FF00FF, 0x0000FFFF };

            c = v;
            c = ((c >> S[0]) & B[0]) + (c & B[0]);
            c = ((c >> S[1]) & B[1]) + (c & B[1]);
            c = ((c >> S[2]) & B[2]) + (c & B[2]);
            c = ((c >> S[3]) & B[3]) + (c & B[3]);
            c = ((c >> S[4]) & B[4]) + (c & B[4]);

            return (c);
        }
        static ulong computeChecksum(ulong source, int nib)
        {
            /*
            Bits 31 to 30 = 2 LSBs of the previous transmitted word, D29* and D30*
            Bits 29 to  6 = Source data bits, d1, d2, ..., d24
            Bits  5 to  0 = Empty parity bits
            */

            /*
            Bits 31 to 30 = 2 LSBs of the previous transmitted word, D29* and D30*
            Bits 29 to  6 = Data bits transmitted by the SV, D1, D2, ..., D24
            Bits  5 to  0 = Computed parity bits, D25, D26, ..., D30
            */

            /*
            1            2           3
            bit    12 3456 7890 1234 5678 9012 3456 7890
            ---    -------------------------------------
            D25    11 1011 0001 1111 0011 0100 1000 0000
            D26    01 1101 1000 1111 1001 1010 0100 0000
            D27    10 1110 1100 0111 1100 1101 0000 0000
            D28    01 0111 0110 0011 1110 0110 1000 0000
            D29    10 1011 1011 0001 1111 0011 0100 0000
            D30    00 1011 0111 1010 1000 1001 1100 0000
            */

            ulong[] bmask = {
        0x3B1F3480UL, 0x1D8F9A40UL, 0x2EC7CD00UL,
        0x1763E680UL, 0x2BB1F340UL, 0x0B7A89C0UL };

            ulong D;
            ulong d = source & 0x3FFFFFC0UL;
            ulong D29 = (source >> 31) & 0x1UL;
            ulong D30 = (source >> 30) & 0x1UL;

            if (nib != 0) // Non-information bearing bits for word 2 and 10
            {
                /*
                Solve bits 23 and 24 to presearve parity check
                with zeros in bits 29 and 30.
                */

                if ((D30 + countBits(bmask[4] & d)) % 2 != 0)
                    d ^= (0x1UL << 6);
                if ((D29 + countBits(bmask[5] & d)) % 2 != 0)
                    d ^= (0x1UL << 7);
            }

            D = d;
            if (D30 != 0)
                D ^= 0x3FFFFFC0UL;

            D |= ((D29 + countBits(bmask[0] & d)) % 2) << 5;
            D |= ((D30 + countBits(bmask[1] & d)) % 2) << 4;
            D |= ((D29 + countBits(bmask[2] & d)) % 2) << 3;
            D |= ((D30 + countBits(bmask[3] & d)) % 2) << 2;
            D |= ((D30 + countBits(bmask[4] & d)) % 2) << 1;
            D |= ((D29 + countBits(bmask[5] & d)) % 2);

            D &= 0x3FFFFFFFUL;
            //D |= (source & 0xC0000000UL); // Add D29* and D30* from source data bits

            return (D);
        }
        static int replaceExpDesignator(char[] str, int len)
        {
            int i, n = 0;

            for (i = 0; i < len; i++)
            {
                if (str[i] == 'D')
                {
                    n++;
                    str[i] = 'E';
                }
            }

            return (n);
        }
        static double subGpsTime(gpstime_t g1, gpstime_t g0)
        {
            double dt;

            dt = g1.sec - g0.sec;
            dt += (double)(g1.week - g0.week) * SECONDS_IN_WEEK;

            return (dt);
        }
        static gpstime_t incGpsTime(gpstime_t g0, double dt)
        {
            gpstime_t g1;

            g1.week = g0.week;
            g1.sec = g0.sec + dt;

            g1.sec = Round(g1.sec * 1000.0) / 1000.0; // Avoid rounding error

            while (g1.sec >= SECONDS_IN_WEEK)
            {
                g1.sec -= SECONDS_IN_WEEK;
                g1.week++;
            }

            while (g1.sec < 0.0)
            {
                g1.sec += SECONDS_IN_WEEK;
                g1.week--;
            }

            return (g1);
        }
        static int readRinexNavAll(ephem_t[][] eph, ref ionoutc_t ionoutc, string fname)
        {

            FileStream fs = new FileStream(fname, FileMode.Open);
            StreamReader sr = new StreamReader(fs, Encoding.ASCII);
            int ieph;

            int sv;
            string str;
            char[] tmp = new char[20];

            datetime_t t = new datetime_t();
            gpstime_t g = new gpstime_t();
            gpstime_t g0 = new gpstime_t();
            double dt;

            int flags = 0x0;

            // Clear valid flag
            for (ieph = 0; ieph < EPHEM_ARRAY_SIZE; ieph++)
                for (sv = 0; sv < MAX_SAT; sv++)
                    eph[ieph][sv].vflg = 0;

            // Read header lines
            while (true)
            {
                str = sr.ReadLine();
                if (str == null)
                    break;

                if (strncmp(str, 60, "END OF HEADER", 13) == true)
                    break;
                else if (strncmp(str, 60, "ION ALPHA", 9) == true)
                {
                    strncpy(tmp, str, 2, 12);
                    tmp[12] = '\0';
                    replaceExpDesignator(tmp, 12);

                    ionoutc.alpha0 = double.Parse(new string(tmp));

                    strncpy(tmp, str, 14, 12);
                    tmp[12] = '\0';
                    replaceExpDesignator(tmp, 12);
                    ionoutc.alpha1 = double.Parse(new string(tmp));

                    strncpy(tmp, str, 26, 12);
                    tmp[12] = '\0';
                    replaceExpDesignator(tmp, 12);
                    ionoutc.alpha2 = double.Parse(new string(tmp));

                    strncpy(tmp, str, 38, 12);
                    tmp[12] = '\0';
                    replaceExpDesignator(tmp, 12);
                    ionoutc.alpha3 = double.Parse(new string(tmp));

                    flags |= 0x1;
                }
                else if (strncmp(str, 60, "ION BETA", 8) == true)
                {
                    strncpy(tmp, str, 2, 12);
                    tmp[12] = '\0';
                    replaceExpDesignator(tmp, 12);
                    ionoutc.beta0 = double.Parse(new string(tmp));

                    strncpy(tmp, str, 14, 12);
                    tmp[12] = '\0';
                    replaceExpDesignator(tmp, 12);
                    ionoutc.beta1 = double.Parse(new string(tmp));

                    strncpy(tmp, str, 26, 12);
                    tmp[12] = '\0';
                    replaceExpDesignator(tmp, 12);
                    ionoutc.beta2 = double.Parse(new string(tmp));

                    strncpy(tmp, str, 38, 12);
                    tmp[12] = '\0';
                    replaceExpDesignator(tmp, 12);
                    ionoutc.beta3 = double.Parse(new string(tmp));

                    flags |= 0x1 << 1;
                }
                else if (strncmp(str, 60, "DELTA-UTC", 9) == true)
                {
                    strncpy(tmp, str, 3, 19);
                    tmp[19] = '\0';
                    replaceExpDesignator(tmp, 19);
                    ionoutc.A0 = double.Parse(new string(tmp));

                    strncpy(tmp, str, 22, 19);
                    tmp[19] = '\0';
                    replaceExpDesignator(tmp, 19);
                    ionoutc.A1 = double.Parse(new string(tmp));

                    strncpy(tmp, str, 41, 9);
                    tmp[9] = '\0';
                    ionoutc.tot = int.Parse(new string(tmp));

                    strncpy(tmp, str, 50, 9);
                    tmp[9] = '\0';
                    ionoutc.wnt = int.Parse(new string(tmp));

                    if (ionoutc.tot % 4096 == 0)
                        flags |= 0x1 << 2;
                }
                else if (strncmp(str, 60, "LEAP SECONDS", 12) == true)
                {
                    strncpy(tmp, str, 0, 6);
                    tmp[6] = '\0';
                    ionoutc.dtls = int.Parse(new string(tmp));

                    flags |= 0x1 << 3;
                }
            }

            ionoutc.vflg = 0;
            if (flags == 0xF) // Read all Iono/UTC lines
                ionoutc.vflg = 1;

            // Read ephemeris blocks
            g0.week = -1;
            ieph = 0;

            while (true)
            {
                str = sr.ReadLine();
                if (str == null)
                    break;

                // PRN
                strncpy(tmp, str, 0, 2);
                tmp[2] = '\0';
                sv = int.Parse(new string(tmp)) - 1;

                // EPOCH
                strncpy(tmp, str, 3, 2);
                tmp[2] = '\0';
                t.y = int.Parse(new string(tmp)) + 2000;

                strncpy(tmp, str, 6, 2);
                tmp[2] = '\0';
                t.m = int.Parse(new string(tmp));

                strncpy(tmp, str, 9, 2);
                tmp[2] = '\0';
                t.d = int.Parse(new string(tmp));

                strncpy(tmp, str, 12, 2);
                tmp[2] = '\0';
                t.hh = int.Parse(new string(tmp));

                strncpy(tmp, str, 15, 2);
                tmp[2] = '\0';
                t.mm = int.Parse(new string(tmp));

                strncpy(tmp, str, 18, 4);
                tmp[4] = '\0';
                t.sec = double.Parse(new string(tmp));

                date2gps(t, ref g);

                if (g0.week == -1)
                    g0 = g;

                // Check current time of clock
                dt = subGpsTime(g, g0);

                if (dt > SECONDS_IN_HOUR)
                {
                    g0 = g;
                    ieph++; // a new set of ephemerides

                    if (ieph >= EPHEM_ARRAY_SIZE)
                        break;
                }

                // Date and time
                eph[ieph][sv].t = t;

                // SV CLK
                eph[ieph][sv].toc = g;

                strncpy(tmp, str, 22, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19); // tmp[15]='E';
                eph[ieph][sv].af0 = double.Parse(new string(tmp));

                strncpy(tmp, str, 41, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].af1 = double.Parse(new string(tmp));

                strncpy(tmp, str, 60, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].af2 = double.Parse(new string(tmp));

                // BROADCAST ORBIT - 1
                str = sr.ReadLine();
                if (str == null)
                    break;

                strncpy(tmp, str, 3, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].iode = (int)double.Parse(new string(tmp));

                strncpy(tmp, str, 22, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].crs = double.Parse(new string(tmp));

                strncpy(tmp, str, 41, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].deltan = double.Parse(new string(tmp));

                strncpy(tmp, str, 60, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].m0 = double.Parse(new string(tmp));

                // BROADCAST ORBIT - 2
                str = sr.ReadLine();
                if (str == null)
                    break;

                strncpy(tmp, str, 3, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].cuc = double.Parse(new string(tmp));

                strncpy(tmp, str, 22, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].ecc = double.Parse(new string(tmp));

                strncpy(tmp, str, 41, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].cus = double.Parse(new string(tmp));

                strncpy(tmp, str, 60, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].sqrta = double.Parse(new string(tmp));

                // BROADCAST ORBIT - 3
                str = sr.ReadLine();
                if (str == null)
                    break;

                strncpy(tmp, str, 3, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].toe.sec = double.Parse(new string(tmp));

                strncpy(tmp, str, 22, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].cic = double.Parse(new string(tmp));

                strncpy(tmp, str, 41, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].omg0 = double.Parse(new string(tmp));

                strncpy(tmp, str, 60, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].cis = double.Parse(new string(tmp));

                // BROADCAST ORBIT - 4
                str = sr.ReadLine();
                if (str == null)
                    break;

                strncpy(tmp, str, 3, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].inc0 = double.Parse(new string(tmp));

                strncpy(tmp, str, 22, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].crc = double.Parse(new string(tmp));

                strncpy(tmp, str, 41, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].aop = double.Parse(new string(tmp));

                strncpy(tmp, str, 60, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].omgdot = double.Parse(new string(tmp));

                // BROADCAST ORBIT - 5
                str = sr.ReadLine();
                if (str == null)
                    break;

                strncpy(tmp, str, 3, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].idot = double.Parse(new string(tmp));

                strncpy(tmp, str, 22, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].codeL2 = (int)double.Parse(new string(tmp));

                strncpy(tmp, str, 41, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].toe.week = (int)double.Parse(new string(tmp));

                // BROADCAST ORBIT - 6
                str = sr.ReadLine();
                if (str == null)
                    break;

                strncpy(tmp, str, 22, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].svhlth = (int)double.Parse(new string(tmp));
                if ((eph[ieph][sv].svhlth > 0) && (eph[ieph][sv].svhlth < 32))

                    eph[ieph][sv].svhlth += 32; // Set MSB to 1

                strncpy(tmp, str, 41, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].tgd = double.Parse(new string(tmp));

                strncpy(tmp, str, 60, 19);
                tmp[19] = '\0';
                replaceExpDesignator(tmp, 19);
                eph[ieph][sv].iodc = (int)double.Parse(new string(tmp));

                // BROADCAST ORBIT - 7
                str = sr.ReadLine();
                if (str == null)
                    break;

                // Set valid flag
                eph[ieph][sv].vflg = 1;

                // Update the working variables
                eph[ieph][sv].A = eph[ieph][sv].sqrta * eph[ieph][sv].sqrta;
                eph[ieph][sv].n = Sqrt(GM_EARTH / (eph[ieph][sv].A * eph[ieph][sv].A * eph[ieph][sv].A)) + eph[ieph][sv].deltan;
                eph[ieph][sv].sq1e2 = Sqrt(1.0 - eph[ieph][sv].ecc * eph[ieph][sv].ecc);
                eph[ieph][sv].omgkdot = eph[ieph][sv].omgdot - OMEGA_EARTH;
            }

            sr.Close();
            fs.Close();

            if (g0.week >= 0)
                ieph += 1; // Number of sets of ephemerides

            return (ieph);
        }
        static private void strncpy(char[] tmp, string str, int v1, int v2)
        {
            for (int i = 0; i < tmp.Length; i++)
                if (i < v2)
                    tmp[i] = str[v1 + i];
                else
                    tmp[i] = '\0';

        }
        static bool strncmp(string v1, int p, string v2, int v3)
        {
            for (int i = 0; i < v2.Length; i++)
                if (v1[p + i] != v2[i]) return false;
            return true;
        }
        static double ionosphericDelay(ionoutc_t ionoutc, gpstime_t g, double[] llh, double[] azel)
        {

            double iono_delay = 0.0;
            double E, phi_u, lam_u, F;

            if (ionoutc.enable == 0)
                return (0.0); // No ionospheric delay

            E = azel[1] / PI;
            phi_u = llh[0] / PI;
            lam_u = llh[1] / PI;

            // Obliquity factor
            F = 1.0 + 16.0 * Pow((0.53 - E), 3.0);

            if (ionoutc.vflg == 0)
                iono_delay = F * 5.0e-9 * SPEED_OF_LIGHT;
            else
            {
                double t, psi, phi_i, lam_i, phi_m, phi_m2, phi_m3;
                double AMP, PER, X, X2, X4;

                // Earth's central angle between the user position and the earth projection of
                // ionospheric intersection point (semi-circles)
                psi = 0.0137 / (E + 0.11) - 0.022;

                // Geodetic latitude of the earth projection of the ionospheric intersection point
                // (semi-circles)
                phi_i = phi_u + psi * Cos(azel[0]);
                if (phi_i > 0.416)
                    phi_i = 0.416;
                else if (phi_i < -0.416)
                    phi_i = -0.416;

                // Geodetic longitude of the earth projection of the ionospheric intersection point
                // (semi-circles)
                lam_i = lam_u + psi * Sin(azel[0]) / Cos(phi_i * PI);

                // Geomagnetic latitude of the earth projection of the ionospheric intersection
                // point (mean ionospheric height assumed 350 km) (semi-circles)
                phi_m = phi_i + 0.064 * Cos((lam_i - 1.617) * PI);
                phi_m2 = phi_m * phi_m;
                phi_m3 = phi_m2 * phi_m;

                AMP = ionoutc.alpha0 + ionoutc.alpha1 * phi_m
                    + ionoutc.alpha2 * phi_m2 + ionoutc.alpha3 * phi_m3;
                if (AMP < 0.0)
                    AMP = 0.0;

                PER = ionoutc.beta0 + ionoutc.beta1 * phi_m
                    + ionoutc.beta2 * phi_m2 + ionoutc.beta3 * phi_m3;
                if (PER < 72000.0)
                    PER = 72000.0;

                // Local time (sec)
                t = SECONDS_IN_DAY / 2.0 * lam_i + g.sec;
                while (t >= SECONDS_IN_DAY)
                    t -= SECONDS_IN_DAY;
                while (t < 0)
                    t += SECONDS_IN_DAY;

                // Phase (radians)
                X = 2.0 * PI * (t - 50400.0) / PER;

                if (Abs(X) < 1.57)
                {
                    X2 = X * X;
                    X4 = X2 * X2;
                    iono_delay = F * (5.0e-9 + AMP * (1.0 - X2 / 2.0 + X4 / 24.0)) * SPEED_OF_LIGHT;
                }
                else
                    iono_delay = F * 5.0e-9 * SPEED_OF_LIGHT;
            }

            return (iono_delay);
        }
        static void computeRange(ref range_t rho, ephem_t eph, ionoutc_t ionoutc, gpstime_t g, double[] xyz)
        {
            double[] pos = new double[3], vel = new double[3], clk = new double[2];
            double[] los = new double[3];
            double tau;
            double range, rate;
            double xrot, yrot;

            double[] llh = new double[3], neu = new double[3];
            double[][] tmat = new double[3][];
            for (int j = 0; j < 3; j++)
                tmat[j] = new double[3];

            // SV position at time of the pseudorange observation.
            satpos(eph, g, pos, vel, clk);

            // Receiver to satellite vector and light-time.
            subVect(los, pos, xyz);
            tau = normVect(los) / SPEED_OF_LIGHT;

            // Extrapolate the satellite position backwards to the transmission time.
            pos[0] -= vel[0] * tau;
            pos[1] -= vel[1] * tau;
            pos[2] -= vel[2] * tau;

            // Earth rotation correction. The change in velocity can be neglected.
            xrot = pos[0] + pos[1] * OMEGA_EARTH * tau;
            yrot = pos[1] - pos[0] * OMEGA_EARTH * tau;
            pos[0] = xrot;
            pos[1] = yrot;

            // New observer to satellite vector and satellite range.
            subVect(los, pos, xyz);
            range = normVect(los);
            rho.d = range;

            // Pseudorange.
            rho.range = range - SPEED_OF_LIGHT * clk[0];

            // Relative velocity of SV and receiver.
            rate = dotProd(vel, los) / range;

            // Pseudorange rate.
            rho.rate = rate; // - SPEED_OF_LIGHT*clk[1];

            // Time of application.
            rho.g = g;

            // Azimuth and elevation angles.
            xyz2llh(xyz, llh);
            ltcmat(llh, tmat);
            ecef2neu(los, tmat, neu);
            neu2azel(rho.azel, neu);

            // Add ionospheric delay
            rho.iono_delay = ionosphericDelay(ionoutc, g, llh, rho.azel);
            rho.range += rho.iono_delay;

            return;
        }
        static void computeCodePhase(ref channel_t chan, range_t rho1, double dt)
        {
            double ms;
            int ims;
            double rhorate;

            // Pseudorange rate.
            rhorate = (rho1.range - chan.rho0.range) / dt;

            // Carrier and code frequency.
            chan.f_carr = -rhorate / LAMBDA_L1;
            chan.f_code = CODE_FREQ + chan.f_carr * CARR_TO_CODE;

            // Initial code phase and data bit counters.
            ms = ((subGpsTime(chan.rho0.g, chan.g0) + 6.0) - chan.rho0.range / SPEED_OF_LIGHT) * 1000.0;

            ims = (int)ms;
            chan.code_phase = (ms - (double)ims) * CA_SEQ_LEN; // in chip

            chan.iword = ims / 600; // 1 word = 30 bits = 600 ms
            ims -= chan.iword * 600;

            chan.ibit = ims / 20; // 1 bit = 20 code = 20 ms
            ims -= chan.ibit * 20;

            chan.icode = ims; // 1 code = 1 ms

            chan.codeCA = chan.ca[(int)chan.code_phase] * 2 - 1;
            chan.dataBit = (int)((chan.dwrd[chan.iword] >> (29 - chan.ibit)) & 0x1UL) * 2 - 1;

            // Save current pseudorange
            chan.rho0 = rho1;

            return;
        }
        static int generateNavMsg(gpstime_t g, ref channel_t chan, int init)
        {
            int iwrd, isbf;
            gpstime_t g0;
            ulong wn, tow;
            ulong sbfwrd;
            ulong prevwrd = 0UL;
            int nib;

            g0.week = g.week;
            g0.sec = (((ulong)(g.sec + 0.5)) / 30UL) * 30.0; // Align with the full frame length = 30 sec
            chan.g0 = g0; // Data bit reference time

            wn = (ulong)(g0.week % 1024);
            tow = ((ulong)g0.sec) / 6UL;

            if (init == 1) // Initialize subframe 5
            {
                prevwrd = 0UL;

                for (iwrd = 0; iwrd < N_DWRD_SBF; iwrd++)
                {
                    sbfwrd = chan.sbf[4][iwrd];

                    // Add TOW-count message into HOW
                    if (iwrd == 1)
                        sbfwrd |= ((tow & 0x1FFFFUL) << 13);

                    // Compute checksum
                    sbfwrd |= (prevwrd << 30) & 0xC0000000UL; // 2 LSBs of the previous transmitted word
                    nib = ((iwrd == 1) || (iwrd == 9)) ? 1 : 0; // Non-information bearing bits for word 2 and 10
                    chan.dwrd[iwrd] = computeChecksum(sbfwrd, nib);

                    prevwrd = chan.dwrd[iwrd];
                }
            }
            else // Save subframe 5
            {
                for (iwrd = 0; iwrd < N_DWRD_SBF; iwrd++)
                {
                    chan.dwrd[iwrd] = chan.dwrd[N_DWRD_SBF * N_SBF + iwrd];

                    prevwrd = chan.dwrd[iwrd];
                }
                /*
                // Sanity check
                if (((chan.dwrd[1])&(0x1FFFFUL<<13)) != ((tow&0x1FFFFUL)<<13))
                {
                Console.WriteLine("\nWARNING: Invalid TOW in subframe 5.\n");
                return(0);
                }
                */
            }

            for (isbf = 0; isbf < N_SBF; isbf++)
            {
                tow++;

                for (iwrd = 0; iwrd < N_DWRD_SBF; iwrd++)
                {
                    sbfwrd = chan.sbf[isbf][iwrd];

                    // Add transmission week number to Subframe 1
                    if ((isbf == 0) && (iwrd == 2))
                        sbfwrd |= (wn & 0x3FFUL) << 20;

                    // Add TOW-count message into HOW
                    if (iwrd == 1)
                        sbfwrd |= ((tow & 0x1FFFFUL) << 13);

                    // Compute checksum
                    sbfwrd |= (prevwrd << 30) & 0xC0000000UL; // 2 LSBs of the previous transmitted word
                    nib = ((iwrd > 0) || (iwrd == 9)) ? 1 : 0; // Non-information bearing bits for word 2 and 10
                    chan.dwrd[(isbf + 1) * N_DWRD_SBF + iwrd] = computeChecksum(sbfwrd, nib);

                    prevwrd = chan.dwrd[(isbf + 1) * N_DWRD_SBF + iwrd];
                }
            }

            return (1);
        }
        static int checkSatVisibility(ephem_t eph, gpstime_t g, double[] xyz, double elvMask, double[] azel)
        {
            double[] llh = new double[3], neu = new double[3];
            double[] pos = new double[3], vel = new double[3], clk = new double[3], los = new double[3];
            double[][] tmat = new double[3][];
            for (int j = 0; j < 3; j++)
                tmat[j] = new double[3];

            if (eph.vflg != 1)
                return (-1); // Invalid

            xyz2llh(xyz, llh);
            ltcmat(llh, tmat);

            satpos(eph, g, pos, vel, clk);
            subVect(los, pos, xyz);
            ecef2neu(los, tmat, neu);
            neu2azel(azel, neu);

            if (azel[1] * R2D > elvMask)
                return (1); // Visible
                            // else
            return (0); // Invisible
        }
        static int allocateChannel(channel_t[] chan, ephem_t[] eph, ionoutc_t ionoutc, gpstime_t grx, double[] xyz, double elvMask)
        {
            int nsat = 0;
            int i, sv;
            double[] azel = new double[2];

            range_t rho = new range_t();
            rho.azel = new double[2];
            double[] refe = { 0.0, 0.0, 0.0 };


            double r_ref, r_xyz;
            double phase_ini;

            for (sv = 0; sv < MAX_SAT; sv++)
            {
                if (checkSatVisibility(eph[sv], grx, xyz, 0.0, azel) == 1)
                {
                    nsat++; // Number of visible satellites

                    if (allocatedSat[sv] == -1) // Visible but not allocated
                    {
                        // Allocated new satellite
                        for (i = 0; i < MAX_CHAN; i++)
                        {
                            if (chan[i].prn == 0)
                            {
                                // Initialize channel
                                chan[i].prn = sv + 1;
                                chan[i].azel[0] = azel[0];
                                chan[i].azel[1] = azel[1];

                                // C/A code generation
                                codegen(chan[i].ca, chan[i].prn);

                                // Generate subframe
                                eph2sbf(eph[sv], ionoutc, chan[i].sbf);

                                // Generate navigation message
                                generateNavMsg(grx, ref chan[i], 1);

                                // Initialize pseudorange
                                computeRange(ref rho, eph[sv], ionoutc, grx, xyz);
                                chan[i].rho0 = rho;

                                // Initialize carrier phase
                                r_xyz = rho.range;

                                computeRange(ref rho, eph[sv], ionoutc, grx, refe);
                                r_ref = rho.range;

                                phase_ini = (2.0 * r_ref - r_xyz) / LAMBDA_L1;
                                chan[i].carr_phase = phase_ini - Floor(phase_ini);

                                // Done.
                                break;
                            }
                        }

                        // Set satellite allocation channel
                        if (i < MAX_CHAN)
                            allocatedSat[sv] = i;
                    }
                }
                else if (allocatedSat[sv] >= 0) // Not visible but allocated
                {
                    // Clear channel
                    chan[allocatedSat[sv]].prn = 0;

                    // Clear satellite allocation flag
                    allocatedSat[sv] = -1;
                }
            }

            return (nsat);
        }
        
        public static double inacc = 0;
        public static double[] Loc { get; set; }
        public static void GPS_SIM(string nav)
        {

            int sv;
            int neph, ieph;
            ephem_t[][] eph = new ephem_t[EPHEM_ARRAY_SIZE][];
            for (int j = 0; j < EPHEM_ARRAY_SIZE; j++)
                eph[j] = new ephem_t[32];
            gpstime_t g0 = new gpstime_t();
            int i;
            channel_t[] chan = new channel_t[MAX_CHAN];
            double elvmask = 0.0; // in degree
            int ip, qp;
            int iTable;
            gpstime_t grx = new gpstime_t();
            double delt;
            int isamp;
            double[] xyz = new double[3];
            double samp_freq = 2.6e6;
            int iq_buff_size;
            int[] gain = new int[MAX_CHAN];
            double path_loss;
            double ant_gain;
            double[] ant_pat = new double[37];
            int ibs; // boresight angle index
            datetime_t t0 = new datetime_t(), tmin = new datetime_t(), tmax = new datetime_t();
            gpstime_t gmin = new gpstime_t(), gmax = new gpstime_t();
            double dt;
            ionoutc_t ionoutc = new ionoutc_t();


            // Default options
            g0.week = -1; // Invalid start time
            ionoutc.enable = 1;
            t0.y = DateTime.UtcNow.Year;
            t0.m = DateTime.UtcNow.Month;
            t0.d = DateTime.UtcNow.Day;
            t0.hh = DateTime.UtcNow.Hour;
            t0.mm = DateTime.UtcNow.Minute;
            t0.sec = (double)DateTime.UtcNow.Second;
            date2gps(t0, ref g0);
            samp_freq = Floor(samp_freq / 10.0);
            iq_buff_size = (int)samp_freq; // samples per 0.1sec
            samp_freq *= 10.0;
            delt = 1.0 / samp_freq;

            ////////////////////////////////////////////////////////////
            // Read ephemeris
            ////////////////////////////////////////////////////////////
            neph = readRinexNavAll(eph, ref ionoutc, nav);  
            for (sv = 0; sv < MAX_SAT; sv++)
            {
                if (eph[0][sv].vflg == 1)
                {
                    gmin = eph[0][sv].toc;
                    tmin = eph[0][sv].t;
                    break;
                }
            }
            for (sv = 0; sv < MAX_SAT; sv++)
            {
                if (eph[neph - 1][sv].vflg == 1)
                {
                    gmax = eph[neph - 1][sv].toc;
                    tmax = eph[neph - 1][sv].t;
                    break;
                }
            }
            if (g0.week >= 0) // Scenario start time has been set.
            {

                gpstime_t gtmp = new gpstime_t();
                datetime_t ttmp = new datetime_t();
                double dsec;

                gtmp.week = g0.week;
                gtmp.sec = (double)(((int)(g0.sec)) / 7200) * 7200.0;

                dsec = subGpsTime(gtmp, gmin);

                // Overwrite the UTC reference week number
                ionoutc.wnt = gtmp.week;
                ionoutc.tot = (int)gtmp.sec;

                // Iono/UTC parameters may no longer valid
                //ionoutc.vflg = FALSE;

                // Overwrite the TOC and TOE to the scenario start time
                for (sv = 0; sv < MAX_SAT; sv++)
                {
                    for (i = 0; i < neph; i++)
                    {
                        if (eph[i][sv].vflg == 1)
                        {
                            gtmp = incGpsTime(eph[i][sv].toc, dsec);
                            gps2date(gtmp, ref ttmp);
                            eph[i][sv].toc = gtmp;
                            eph[i][sv].t = ttmp;

                            gtmp = incGpsTime(eph[i][sv].toe, dsec);
                            eph[i][sv].toe = gtmp;
                        }
                    }
                }

            }
            else
            {
                g0 = gmin;
                t0 = tmin;
            }
            // Select the current set of ephemerides
            ieph = -1;
            for (i = 0; i < neph; i++)
            {
                for (sv = 0; sv < MAX_SAT; sv++)
                {
                    if (eph[i][sv].vflg == 1)
                    {
                        dt = subGpsTime(g0, eph[i][sv].toc);
                        if (dt >= -SECONDS_IN_HOUR && dt < SECONDS_IN_HOUR)
                        {
                            ieph = i;
                            break;
                        }
                    }
                }

                if (ieph >= 0) // ieph has been set
                    break;
            }
            if (ieph == -1)
            {
                return;
            }

            ////////////////////////////////////////////////////////////
            // Allocate I/Q buffer
            ////////////////////////////////////////////////////////////
            short[] iq_buff = new short[2 * iq_buff_size];


            ////////////////////////////////////////////////////////////
            // Initialize channels
            ////////////////////////////////////////////////////////////
            // Clear all channels
            for (i = 0; i < MAX_CHAN; i++)
            {
                chan[i].prn = 0;
                chan[i].azel = new double[2];
                chan[i].ca = new int[CA_SEQ_LEN];
                chan[i].sbf = new ulong[5][];
                for (int j = 0; j < 5; j++)
                    chan[i].sbf[j] = new ulong[10];
                chan[i].dwrd = new ulong[60];
            }
            // Clear satellite allocation flag
            for (sv = 0; sv < MAX_SAT; sv++)
                allocatedSat[sv] = -1;
            // Initial reception time
            grx = incGpsTime(g0, 0.0);
            llh2xyz(Loc, xyz);
            // Allocate visible satellites
            allocateChannel(chan, eph[ieph], ionoutc, grx, xyz, elvmask);

            ////////////////////////////////////////////////////////////
            // Receiver antenna gain pattern
            ////////////////////////////////////////////////////////////
            for (i = 0; i < 37; i++)
                ant_pat[i] = Pow(10.0, -ant_pat_db[i] / 20.0);


            ////////////////////////////////////////////////////////////
            // Generate baseband signals
            ////////////////////////////////////////////////////////////
            // Update receiver time
            grx = incGpsTime(grx, 0.1);
            Stopwatch start = new Stopwatch();
            start.Start();

            while (true)
            {
                for (i = 0; i < MAX_CHAN; i++)
                {
                    if (chan[i].prn > 0)
                    {
                        // Refresh code phase and data bit counters
                        range_t rho = new range_t();
                        rho.azel = new double[2];
                        sv = chan[i].prn - 1;


                        computeRange(ref rho, eph[ieph][sv], ionoutc, grx, xyz);


                        chan[i].azel[0] = rho.azel[0];
                        chan[i].azel[1] = rho.azel[1];

                        // Update code phase and data bit counters
                        computeCodePhase(ref chan[i], rho, 0.1);
                        chan[i].carr_phasestep = (int)(512 * 65536.0 * chan[i].f_carr * delt);

                        // Path loss
                        path_loss = 20200000.0 / rho.d;

                        // Receiver antenna gain
                        ibs = (int)((90.0 - rho.azel[1] * R2D) / 5.0); // covert elevation to boresight
                        ant_gain = ant_pat[ibs];

                        // Signal gain
                        gain[i] = (int)(path_loss * ant_gain * 128.0); // scaled by 2^7
                    }
                }
                for (isamp = 0; isamp < iq_buff_size; isamp++)
                {
                    int i_acc = 0;
                    int q_acc = 0;

                    for (i = 0; i < MAX_CHAN; i++)
                    {
                        if (chan[i].prn > 0)
                        {
                            iTable = ((int)chan[i].carr_phase >> 16) & 511;

                            ip = chan[i].dataBit * chan[i].codeCA * cosTable512[iTable] * gain[i];
                            qp = chan[i].dataBit * chan[i].codeCA * sinTable512[iTable] * gain[i];

                            // Accumulate for all visible satellites
                            i_acc += ip;
                            q_acc += qp;

                            // Update code phase
                            chan[i].code_phase += chan[i].f_code * delt;

                            if (chan[i].code_phase >= CA_SEQ_LEN)
                            {
                                chan[i].code_phase -= CA_SEQ_LEN;

                                chan[i].icode++;

                                if (chan[i].icode >= 20) // 20 C/A codes = 1 navigation data bit
                                {
                                    chan[i].icode = 0;
                                    chan[i].ibit++;

                                    if (chan[i].ibit >= 30) // 30 navigation data bits = 1 word
                                    {
                                        chan[i].ibit = 0;
                                        chan[i].iword++;
                                        /*
                                        if (chan[i].iword>=N_DWRD)
                                        Console.WriteLine("\nWARNING: Subframe word buffer overflow.\n");
                                        */
                                    }

                                    // Set new navigation data bit
                                    chan[i].dataBit = (int)((chan[i].dwrd[chan[i].iword] >> (29 - chan[i].ibit)) & 0x1UL) * 2 - 1;
                                }
                            }

                            // Set currnt code chip
                            chan[i].codeCA = chan[i].ca[(int)chan[i].code_phase] * 2 - 1;

                            // Update carrier phase
                            chan[i].carr_phase += chan[i].carr_phasestep;
                        }
                    }

                    // Scaled by 2^7
                    i_acc = (i_acc + 64) >> 7;
                    q_acc = (q_acc + 64) >> 7;

                    // Store I/Q samples into buffer

                    iq_buff[isamp * 2] = (short)i_acc;
                    iq_buff[isamp * 2 + 1] = (short)q_acc;

                }
                byte[] datap = new byte[iq_buff_size * 2];
                for (int l = 0; l < iq_buff_size * 2; l++)
                {
                    datap[l] = (byte)(iq_buff[l] >> 4);
                }
                TCP.TCPsend(datap);

                ////////////////////////////////////////////////////////////
                // Update navigation message and channel allocation
                ////////////////////////////////////////////////////////////

                // Update navigation message
                for (i = 0; i < MAX_CHAN; i++)
                {
                    if (chan[i].prn > 0)
                        generateNavMsg(grx, ref chan[i], 0);
                }

                // Refresh ephemeris and subframes
                // Quick and dirty fix. Need more elegant way.
                for (sv = 0; sv < MAX_SAT; sv++)
                {
                    if (eph[ieph + 1][sv].vflg == 1)
                    {
                        dt = subGpsTime(eph[ieph + 1][sv].toc, grx);
                        if (dt < SECONDS_IN_HOUR)
                        {
                            ieph++;

                            for (i = 0; i < MAX_CHAN; i++)
                            {
                                // Generate new subframes if allocated
                                if (chan[i].prn != 0)
                                    eph2sbf(eph[ieph][chan[i].prn - 1], ionoutc, chan[i].sbf);
                            }
                        }

                        break;
                    }
                }

                // Update channel allocation
                llh2xyz(Loc, xyz);

                //if (subGpsTime(grx, g0) - start.ElapsedMilliseconds / 1000 > 0.1)
                //    Thread.Sleep(10);

                // Update receiver time
                grx = incGpsTime(grx, 0.1);
            }
        }

    }
}



