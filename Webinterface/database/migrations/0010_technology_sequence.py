# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('database', '0009_auto_20150628_1305'),
    ]

    operations = [
        migrations.AddField(
            model_name='technology',
            name='sequence',
            field=models.CharField(default=b'000000', max_length=500),
        ),
    ]
