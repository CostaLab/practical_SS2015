# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('database', '0011_technology_differgraph'),
    ]

    operations = [
        migrations.AlterField(
            model_name='technology',
            name='snp_file',
            field=models.FileField(null=True, upload_to=b'documents', blank=True),
        ),
    ]
