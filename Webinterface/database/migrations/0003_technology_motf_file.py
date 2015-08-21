# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('database', '0002_auto_20150627_1041'),
    ]

    operations = [
        migrations.AddField(
            model_name='technology',
            name='motf_file',
            field=models.FileField(default=b'documents/default.txt', upload_to=b'documents'),
        ),
    ]
