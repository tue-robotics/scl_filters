#include <scl/filters/DFirstOrderLowpass.hpp>
#include <scl/filters/DSecondOrderLowpass.hpp>
#include <scl/filters/DLeadLag.hpp>
#include <scl/filters/DPD.hpp>
#include <scl/filters/DPID.hpp>
#include <scl/filters/DSkewedNotch.hpp>
#include <scl/filters/DWeakIntegrator.hpp>

#include <cmath>
#include <string>
#include <iostream>
#include <vector>

// ----------------------------------------------------------------------------------------------------

double runFilter(DFILTERS::Filter& filter)
{
    for(unsigned int i = 0; i < 1000; ++i)
    {
        // Arbitrarily chosen signal
        double input = 6 * sin(0.001 * i) + 2 * sin(0.05 * i) - 5 * sin(0.3 * i);
        filter.update(input);
    }

    return filter.getOutput();
}

// ----------------------------------------------------------------------------------------------------

bool testFilter(const std::string& name, DFILTERS::Filter& filter, double expected_outcome)
{
    double output = runFilter(filter);

    bool ok = std::abs(filter.getOutput() - expected_outcome) < 1e-5;
    if (ok)
        std::cout << "[OK]   ";
    else
        std::cout << "[ERROR]";


    std::cout << " " << name;

    if (!ok)
    {
        std::cout.precision(10);
        std::cout << "    (outcome = " << filter.getOutput() << ", expected = " << expected_outcome << ")";
    }
    std::cout << std::endl;

    return ok;
}

// ----------------------------------------------------------------------------------------------------

struct FilterInfo
{
    FilterInfo(const std::string& name_, DFILTERS::Filter* filter_) : name(name_), filter(filter_) {}
    std::string name;
    DFILTERS::Filter* filter;
};

// ----------------------------------------------------------------------------------------------------

bool testFilters(std::vector<FilterInfo>& filters, const double* outcomes, int num_methods)
{
    bool ok = true;

    int num_filters = filters.size() / num_methods;

    for(unsigned int i = 0; i < filters.size(); ++i)
    {
        if (i % num_filters == 0)
        {
            std::cout << std::endl;
            std::cout << "--------------------------------------------" << std::endl;
            std::cout << " Discretization method " << (i / num_filters + 1) << std::endl;
            std::cout << "--------------------------------------------" << std::endl;
            std::cout << std::endl;
        }

        FilterInfo& info = filters[i];
        ok &= testFilter(info.name, *info.filter, outcomes[i]);
    }

    return ok;
}

// ----------------------------------------------------------------------------------------------------

void printFilterOutcomes(std::vector<FilterInfo>& filters)
{
    std::cout.precision(10);
    std::cout << "double OUTCOMES[] = {";

    for(unsigned int i = 0; i < filters.size(); ++i)
    {
        FilterInfo& info = filters[i];
        double output = runFilter(*info.filter);
        if (output != output)
            output = 0;

        if (i % 6 == 0)
            std::cout << std::endl << "    ";

        std::cout << output;

        if (i + 1 != filters.size())
            std::cout << ", ";
    }
    std::cout << std::endl << "    };" << std::endl;
}

// ----------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
    std::vector<FilterInfo> filters;

    double dt = 0.001;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Expected outcomes

    // Use printFilterOutcomes below to generate these values
    double OUTCOMES[] = {
        7.443527779, 2.821451805, 111.6014839, 31165.09039, 1229.919097, 13.42303696,
        9.683982941, 7.648917597, 2.462471592, 114.258704, 33964.90763, 1323.218834,
        13.97529244, 9.682254833, 7.506530392, 2.638142232, 114.0757267, 32512.97785,
        1274.834927, 13.7168207, 9.683118887, 7.572172065, 2.637243742, 114.0451243,
        32512.97784, 1274.834465, 13.67293935, 9.683118888, 6.916870485, 9.167897446,
        134.9707448, 34140.50203, 1329.071981, 9.167897446, 9.682254833, 6.916870485,
        9.167897446, 150.7131111, 32489.6856, 9.167897446, 9.167897446, 9.682206359
        };

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Initialize filters

    for(int method = 1; method <= 6; ++method)
    {
        filters.push_back(FilterInfo("DFirstOrderLowpass", new DFILTERS::DFirstOrderLowpass(100, dt, method)));
        filters.push_back(FilterInfo("DSecondOrderLowpass", new DFILTERS::DSecondOrderLowpass(20, 0.7, dt, method)));
        filters.push_back(FilterInfo("DLeadLag", new DFILTERS::DLeadLag(1.6, 60, dt, method)));
        filters.push_back(FilterInfo("DPD", new DFILTERS::DPD(1.6, 60, dt, method)));
        filters.push_back(FilterInfo("DPID", new DFILTERS::DPID(20, 2, 3, dt, method)));
        filters.push_back(FilterInfo("DSkewedNotch", new DFILTERS::DSkewedNotch(60, 1.5, 80, 1.0, dt, method)));
        filters.push_back(FilterInfo("DWeakIntegrator", new DFILTERS::DWeakIntegrator(0.03, dt, method)));
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Use this to generate the outcome array above

//    printFilterOutcomes(filters);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Test the filters

    if (!testFilters(filters, OUTCOMES, 6))
    {
        std::cout << std::endl << "[ERROR]" << std::endl;
        return 1;
    }

    std::cout << std::endl << "[OK]" << std::endl;

    return 0;
}
