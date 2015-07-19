# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
        ('database', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='gene',
            name='refFile',
            field=models.FileField(default=b'documents/default.txt', upload_to=b'documents'),
        ),
        migrations.AddField(
            model_name='gene',
            name='snp_file',
            field=models.FileField(default=b'documents/default.txt', upload_to=b'documents'),
        ),
        migrations.AddField(
            model_name='technology',
            name='barGraph',
            field=models.ImageField(default=b'graphics/default.png', upload_to=b'graphics'),
        ),
        migrations.AddField(
            model_name='technology',
            name='sequence',
            field=models.CharField(default=b'000000', max_length=500),
        ),
    ]
