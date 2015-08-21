# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('database', '0010_technology_sequence'),
    ]

    operations = [
        migrations.AddField(
            model_name='technology',
            name='differGraph',
            field=models.ImageField(default=b'-----', upload_to=b'graphics'),
        ),
    ]
