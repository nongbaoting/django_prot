# Generated by Django 2.2.5 on 2022-03-31 20:59

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0010_cdd'),
    ]

    operations = [
        migrations.CreateModel(
            name='NrCD',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('protin_id', models.CharField(db_index=True, max_length=200)),
                ('length', models.IntegerField(db_index=True)),
                ('desc', models.CharField(db_index=True, max_length=100)),
                ('pssm', models.CharField(db_index=True, max_length=100)),
                ('start', models.IntegerField(db_index=True)),
                ('end', models.IntegerField(db_index=True)),
                ('evalue', models.CharField(db_index=True, max_length=100)),
                ('biscore', models.CharField(db_index=True, max_length=100)),
                ('cdd_name', models.CharField(db_index=True, max_length=1000)),
                ('cdd_nameCat', models.TextField()),
                ('cdd_noteCat', models.TextField()),
                ('cdd_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='protein.CDD')),
            ],
        ),
    ]
