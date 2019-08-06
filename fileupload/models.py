from __future__ import unicode_literals

from django.db import models

class Photo2(models.Model):
    title = models.CharField(max_length=255, blank=True)
    file = models.FileField(upload_to='temp/%Y%m%d%H%M')
    uploaded_at = models.DateTimeField(auto_now_add=True)
    job_id = models.CharField(max_length=255)

