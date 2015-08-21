# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('database', '0007_auto_20150627_1537'),
    ]

    operations = [
        migrations.AddField(
            model_name='technology',
            name='commonGraph',
            field=models.ImageField(default=b'------', upload_to=b'graphics'),
        ),
    ]
