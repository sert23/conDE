from django.conf.urls import url

from . import views

urlpatterns = [
    # url(r'^$', views.new_upload, name='multi_new'),
    url(r'^ajax_graph$', views.test_graph_ajax, name='ajax_graph'),
    url(r'^[A-za-z0-9]+', views.DEresult.as_view(), name="DEresult"),




            ]
