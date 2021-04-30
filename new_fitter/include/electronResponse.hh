#ifndef _electronResponse_h
#define _electronResponse_h

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
class electronResponse
{
    public:
        electronResponse();
        ~electronResponse();

        static double getElecNonl(double Etrue);
        
        static void loadSimElecNonl();
        static void loadMinSimElecNonl();
        static void loadElecResol();
        static void EmpiricalFit();

        static void Plot();

        static void FitPlot();

        static void FuncConstruct();

    private:
        static const int m_nSimData = 809;
        static double m_SimEtrue[m_nSimData];
        static double m_SimNonl[m_nSimData];

        static double m_scale;

        static bool m_loadSimFile;
        static bool m_doFit;
        static bool m_loadResol;

    private:
        static TGraph* gSimData;
        static TF1* fElecNonl;

    private:
        static double m_p0;
        static double m_p1;
        static double m_p2;
        static double m_p3;

        static double m_ra;
        static double m_rb;
        static double m_rc;

    public:
        static double getp0()       {return m_p0;}
        static void setp0(double p0){m_p0 = p0;}
        static double getp1()       {return m_p1;}
        static void setp1(double p1){m_p1 = p1;}
        static double getp2()       {return m_p2;}
        static void setp2(double p2){m_p2 = p2;}
        static double getp3()       {return m_p3;}
        static void setp3(double p3){m_p3 = p3;}
        static bool getLoadResol()  {return m_loadResol;}
        static double getra()       {return m_ra;}
        static void setra(double ra){m_ra = ra;}
        static double getrb()       {return m_rb;}
        static void setrb(double rb){m_rb = rb;}
        static double getrc()       {return m_rc;}
        static void setrc(double rc){m_rc = rc;}

        static void SetParameters();
        static double calcElecNonl(double E)   {return fElecNonl->Eval(E);}

        static TGraphErrors* gMinElecNonl;
        static TGraphErrors* gElecResol;
        static TF1* fElecResol;

};

extern double gElecResolFunc(double* x, double* p);

#endif
