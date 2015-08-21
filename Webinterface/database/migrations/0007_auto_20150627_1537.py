# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('database', '0006_auto_20150627_1345'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='gene',
            name='snp_file',
        ),
        migrations.AddField(
            model_name='technology',
            name='snp_file',
            field=models.FileField(default=b'------', upload_to=b'documents'),
        ),
    ]
