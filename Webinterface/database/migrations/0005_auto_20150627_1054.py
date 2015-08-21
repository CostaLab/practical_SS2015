# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('database', '0004_auto_20150627_1046'),
    ]

    operations = [
        migrations.AlterField(
            model_name='gene',
            name='refFile',
            field=models.FileField(default=b'------', upload_to=b'documents'),
        ),
        migrations.AlterField(
            model_name='gene',
            name='snp_file',
            field=models.FileField(default=b'------', upload_to=b'documents'),
        ),
        migrations.AlterField(
            model_name='technology',
            name='barGraph',
            field=models.ImageField(default=b'------', upload_to=b'graphics'),
        ),
        migrations.AlterField(
            model_name='technology',
            name='motif_file',
            field=models.FileField(default=b'------', upload_to=b'documents'),
        ),
    ]
