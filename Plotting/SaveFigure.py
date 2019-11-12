import os
import datetime

from bokeh import models
from bokeh.io import export_svgs
from holoext.utils import DEFAULT_N, get_cmap, flatten, tidy_fn

def _get_figures_core(objs):
    if isinstance(objs, list):
        objs = [_get_figures_core(plot) for plot in objs]
    elif isinstance(objs, (models.Column, models.Row)):
        objs = [_get_figures_core(child) for child in objs.children
                if not isinstance(child, (models.ToolbarBox,
                                          models.WidgetBox))]
    return objs

def _get_figures(objs):
    try:
        return list(flatten(_get_figures_core(objs)))
    except TypeError:
        return [_get_figures_core(objs)]

def save_to_svg(hv_obj, save=None):
    if save is None:
        save =datetime.date.today().strftime("%Y%m%d")
    bokeh_obj = hv.renderer('bokeh').get_plot(hv_obj).state
    figures = _get_figures(bokeh_obj)
    for i, figure in enumerate(figures):
        figure.output_backend = 'svg'

        if len(figures) != 1:
            if not os.path.exists(save):
                os.mkdir(save)
            tidied_title = tidy_fn(figure.title.text)
            save_fp = os.path.join(
                save, '{0}_{1}'.format(tidied_title, i))
        else:
            save_fp = save

        if not save_fp.endswith('svg'):
            save_fp = '{0}.{1}'.format(save_fp, 'svg')

        export_svgs(figure, save_fp)