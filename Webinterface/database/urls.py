from django.conf.urls import patterns, include, url
from . import views
from django.conf import settings
from django.conf.urls.static import static

from . import views

urlpatterns = [
    url(r'^$', views.gene_list),
    url(r'^gene/(?P<gene_id>[0-9]+)/$', views.gene_detail),
    url(r'^technology_list/(?P<gene_id>[0-9]+)/$', views.technology_list),
    url(r'^technology/(?P<tech_id>[0-9]+)/$', views.technology_detail),
    url(r'^search/$', views.search),
    url(r'^common_search/$', views.common_search)
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)