from django.conf.urls import url

from . import views

urlpatterns = [
    # url(r'^$', views.new_upload, name='multi_new'),
    url(r'^heatmap/recalc$', views.heatmap_recalculate, name='hm_recalc'),
    url(r'^heatmap/[A-za-z0-9]+', views.Heatmap.as_view(), name="heatmap"),


            ]
