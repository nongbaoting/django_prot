# Generated by Django 2.2.5 on 2022-04-12 22:33

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0052_delete_scopecla'),
    ]

    operations = [
        migrations.CreateModel(
            name='ScopeCla',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('domain', models.CharField(db_index=True, max_length=100)),
                ('TP', models.CharField(db_index=True, max_length=100)),
                ('CL', models.CharField(db_index=True, max_length=100)),
                ('CF', models.CharField(db_index=True, max_length=100)),
                ('SF', models.CharField(db_index=True, max_length=100)),
                ('FA', models.CharField(db_index=True, max_length=100)),
                ('protein_type', models.CharField(db_index=True, max_length=200)),
                ('class_fold', models.CharField(db_index=True, max_length=200)),
                ('superfamily', models.CharField(db_index=True, max_length=200)),
                ('family', models.CharField(db_index=True, max_length=200)),
            ],
        ),
    ]
