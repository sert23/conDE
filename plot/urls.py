from django.conf.urls import url

from . import views

urlpatterns = [
    # url(r'^$', views.new_upload, name='multi_new'),

    url(r'^random/[A-za-z0-9]+$', views.random_redirect),
    url(r'^random$', views.heatmap_recalculate, name='random_url'),
    url(r'^heatmap/recalc$', views.heatmap_recalculate, name='hm_recalc'),
    url(r'^heatmap/[A-za-z0-9]+', views.Heatmap.as_view(), name="heatmap"),
    url(r'^volcano/recalc$', views.volcano_recalculate, name='volc_recalc'),
    url(r'^volcano/[A-za-z0-9]+', views.Volcano.as_view(), name="volcano"),


            ]
