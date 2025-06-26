#ifndef CUSTOMGRAPH_H
#define CUSTOMGRAPH_H

#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <string>

namespace KLOE
{

  template <typename T>
  class CustomGraph
  {
  private:
    TGraph *
        graph;
    std::vector<T>
        x_values,
        y_values;

  public:
    CustomGraph(const std::vector<T> &x, const std::vector<T> &y)
    {
      if (x.size() != y.size())
      {
        std::cerr << "Error: x and y vectors must have the same size." << std::endl;
        graph = nullptr;
        return;
      }
      x_values = x;
      y_values = y;
      graph = new TGraph(x.size(), x.data(), y.data());
    }

    ~CustomGraph()
    {
      delete graph;
    }

    void SetStyle(int markerStyle = 20, int lineColor = kBlue, int markerColor = kRed, int lineWidth = 2);

    void Draw(const char *drawOption = "AP");

    void SaveAsPNG(const std::string &filename);
  };
}

// Example usage:
// std::vector<double> x = {1, 2, 3, 4, 5};
// std::vector<double> y = {1, 4, 9, 16, 25};
// CustomGraph<double> graph(x, y);
// graph.SetStyle(21, kGreen, kBlack, 3);
// graph.Draw();
// graph.SaveAsPNG("graph_output");

#endif // !CUSTOMGRAPH_H
