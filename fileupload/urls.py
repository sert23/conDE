from django.conf.urls import url

from . import views

urlpatterns = [
    # url(r'^$', views.new_upload, name='multi_new'),

    # url(r'^random/[A-za-z0-9]+$', views.random_redirect),
    # url(r'^random$', views.heatmap_recalculate, name='random_url'),
    url(r'^$', views.Upload.as_view(), name='uploadf'),
    url(r'groups/[A-za-z0-9]+', views.ChooseGroups.as_view(), name='upload_groups'),
    url(r'groups', views.ChooseGroups.as_view(), name='groups'),
            ]
