from django.conf.urls import url

from . import views

urlpatterns = [
    # url(r'^$', views.new_upload, name='multi_new'),
    url(r'^ajax_graph$', views.get_plot_content, name='ajax_graph'),
    url(r'^ajax_recalc$', views.ajax_recalculate, name='ajax_recalc'),
    url(r'^[A-za-z0-9]+', views.DEresult.as_view(), name="DEresult"),
    url(r'^', views.DEresult.as_view(), name="result"),




            ]
