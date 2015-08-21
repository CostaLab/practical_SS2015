# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('database', '0005_auto_20150627_1054'),
    ]

    operations = [
        migrations.AlterField(
            model_name='technology',
            name='gene',
            field=models.ForeignKey(related_name='technologies', to='database.Gene'),
        ),
    ]
