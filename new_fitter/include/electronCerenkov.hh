#ifndef _ELECTRONCERENKOV_H
#define _ELECTRONCERENKOV_H

#include <vector>
#include <TGraph.h>

using namespace std;

class electronCerenkov
{
    public:
        electronCerenkov();
        ~electronCerenkov();
    
    public:
        static void setkC(double val)           {m_kC = val;}
        static double getkC()                   {return m_kC;}
        static void setEnergyScale(double val)  {m_energyScale = val;}
        static double getEnergyScale()          {return m_energyScale;}

        static void LoadCerenkov();
        
        static double getCerenkovPE(double E);

        static double getCerPE(double E);

        static void Plot();

    private:
        static double m_kC;
        static double m_energyScale;

        static bool m_LoadCerenkov;

        static std::vector<double> m_Etrue;
        static std::vector<double> m_Cerenkov;

        static const int m_nData = 900;
        static double m_E[m_nData];
        static double m_nonl[m_nData];

        static TGraph* gNPE_elec;
};

#endif
