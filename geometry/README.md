# Геометрия скруббера

Только газовая полость: люки, штуцера, проушины, фланцы и изоляция не моделируются.

- `params.py` — все размеры (мм), снятые с чертежа; допущения помечены комментариями.
- `build_geometry.py` — строит модель (CadQuery) и пишет в `out/`:
  - `venturi_fluid.step` — объём газа, мм (открыть в FreeCAD для проверки);
  - `nozzle_solid.step` — форсунка отдельным телом;
  - `inlet.stl`, `outlet.stl`, `wall.stl`, `nozzle.stl` — поверхности в метрах, по одному патчу на файл, для snappyHexMesh.
- `check_drawing.py` — проверочный чертёж `out/check_drawing.png`.

```bash
pip install cadquery matplotlib
python3 build_geometry.py && python3 check_drawing.py
```
