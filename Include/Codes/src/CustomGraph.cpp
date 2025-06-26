#include <CustomGraph.h>

namespace KLOE
{
  template <typename T>
  void CustomGraph<T>::SetStyle(int markerStyle, int lineColor, int markerColor, int lineWidth)
  {
    if (!graph)
      return;
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerColor(markerColor);
    graph->SetLineColor(lineColor);
    graph->SetLineWidth(lineWidth);
  }

  template <typename T>
  void CustomGraph<T>::Draw(const char *drawOption)
  {
    if (!graph)
      return;
    TCanvas *canvas = new TCanvas("canvas", "Graph Canvas", 800, 600);
    graph->Draw(drawOption);
    canvas->Update();
  }

  template <typename T>
  void CustomGraph<T>::SaveAsPNG(const std::string &filename)
  {
    if (!graph)
      return;
    TCanvas *canvas = new TCanvas("canvas", "Graph Canvas", 800, 600);
    graph->Draw("APL");
    canvas->SaveAs((filename + ".png").c_str());
    delete canvas;
  }
}